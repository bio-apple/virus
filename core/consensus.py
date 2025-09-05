import pandas as pd
import matplotlib.pyplot as plt
import subprocess, argparse
import os, random
import numpy as np
import threading
from concurrent.futures import ThreadPoolExecutor

docker = "virus:latest"

def get_blast_description(key, blast_db_dir, blast_db_name, docker, outdir):
    cmd = f'docker run --rm -v {outdir}:/outdir -v {blast_db_dir}:/ref {docker} sh -c \''
    blast_cmd = cmd + f"export PATH=/opt/conda/envs/kraken2/bin/:$PATH && blastdbcmd -db /ref/{blast_db_name} -entry {key} -outfmt \"%a|%t\"\'"
    process = subprocess.Popen(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = process.communicate()

    description = ""
    if process.returncode == 0 and stdout.strip():
        all_descriptions = [line.strip() for line in stdout.strip().split('\n') if line.strip()]
        if all_descriptions:
            description = random.choice(all_descriptions)
    return description


def process_accession(accession, pe1, pe2, outdir, blast_db_dir, blast_db_name, read_length, plot_threshold,fa_threshold, docker, descriptions_lock, descriptions):

    print(f"Starting processing for accession: {accession}")

    # Get description for the accession
    desc = get_blast_description(accession, blast_db_dir, blast_db_name, docker, outdir)
    with descriptions_lock:
        descriptions[accession] = desc

    # Step 1: Get FASTA from BLAST DB
    fasta_cmd = f'docker run --rm -v {outdir}:/outdir -v {blast_db_dir}:/ref {docker} sh -c \'export PATH=/opt/conda/envs/kraken2/bin/:$PATH && blastdbcmd -db /ref/{blast_db_name} -entry {accession} > /outdir/{accession}.fa\''
    subprocess.check_call(fasta_cmd, shell=True)

    # Rename and format the fasta file
    command = ['awk', f'/^>/ {{print ">{accession}"}} !/^>/ {{print}}', f'{outdir}/{accession}.fa']
    with open(f'{outdir}/{accession}.fasta', 'w') as outfile:
        subprocess.run(command, check=True, stdout=outfile)

    # Step 2: Map reads and generate coverage file
    print(f"Mapping reads for {accession}...")
    if pe2 is not None:
        mapping_cmd = (f'docker run --rm '
                       f'-v {pe1}:/raw_data/{os.path.basename(pe1)} '
                       f'-v {pe2}:/raw_data/{os.path.basename(pe2)} '
                       f'-v {outdir}:/outdir/ {docker} '
                       f'sh -c \'export PATH=/opt/conda/bin/:$PATH && cd /outdir/ && '
                       f'bowtie2-build {accession}.fasta {accession}.fasta && '
                       f'bowtie2 --threads 64 -x {accession}.fasta -1 /raw_data/{os.path.basename(pe1)} -2 /raw_data/{os.path.basename(pe2)} |'
                       f'samtools view -bh |samtools sort > /outdir/{accession}.bam && samtools index /outdir/{accession}.bam && '
                       f'pileup.sh in=/outdir/{accession}.bam out=/outdir/{accession}.cov.txt && rm -rf /outdir/{accession}.sam\'')
    else:
        mapping_cmd = (f'docker run --rm '
                       f'-v {pe1}:/raw_data/{os.path.basename(pe1)} '
                       f'-v {outdir}:/outdir/ {docker} '
                       f'sh -c \'export PATH=/opt/conda/bin/:$PATH && cd /outdir/ && '
                       f'bowtie2-build {accession}.fasta {accession}.fasta && '
                       f'bowtie2 --threads 64 -x {accession}.fasta -U /raw_data/{os.path.basename(pe1)} |'
                       f'samtools view -bh |samtools sort > /outdir/{accession}.bam && samtools index /outdir/{accession}.bam && '
                       f'pileup.sh in=/outdir/{accession}.bam out=/outdir/{accession}.cov.txt && rm -rf /outdir/{accession}.sam\'')
    print(mapping_cmd)
    subprocess.check_call(mapping_cmd, shell=True)

    # Step 3: Check coverage and perform further analysis if thresholds are met
    with open(f"{outdir}/{accession}.cov.txt", "r") as f_cov:
        lines_to_write = []
        is_plot_worthy = False
        is_fa_worthy = False

        for line in f_cov:
            line = line.strip()
            if not line.startswith("#"):
                array = line.split("\t")
                if float(array[5]) > max(int(read_length) * 3, 500) or float(array[1]) >= 10 or float(array[9]) >= 10:
                    lines_to_write.append(
                        f"{descriptions.get(array[0], 'N/A')}\t{array[1]}\t{array[2]}\t{array[3]}\t{array[4]}\t{array[5]}\t{array[6]}\t{array[7]}\t{array[8]}\t{array[9]}\n")
                    if float(array[4]) >= plot_threshold:
                        is_plot_worthy = True
                    if float(array[4]) >= fa_threshold:
                        is_fa_worthy = True

        # Perform variant calling and consensus for the accession if fa_threshold is met
        if is_fa_worthy:
            print(f"Variant calling and consensus for {accession}...")
            # Run samtools depth to create mask bed file
            depth_cmd = f'docker run --rm -v {outdir}:/outdir {docker} sh -c \'export PATH=/opt/conda/bin:$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa /outdir/{accession}.bam > /outdir/{accession}.depth.txt\''
            subprocess.check_call(depth_cmd, shell=True)

            with open(f"{outdir}/{accession}.depth.txt", "r") as infile, open(f"{outdir}/{accession}.mask.bed", "w") as outbed:
                for line in infile:
                    chrom, pos, depth = line.strip().split("\t")
                    if int(depth) < 10:
                        outbed.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")

            # Run bcftools mpileup and call
            vcf_cmd = f'docker run --rm -v {outdir}:/outdir {docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bcftools mpileup -Ou -f /outdir/{accession}.fasta /outdir/{accession}.bam | bcftools call --ploidy 1 -mv -Oz -o /outdir/{accession}.vcf.gz\''
            subprocess.check_call(vcf_cmd, shell=True)

            # Run bcftools consensus
            consensus_cmd = f'docker run --rm -v {outdir}:/outdir {docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bcftools index /outdir/{accession}.vcf.gz && cat /outdir/{accession}.fasta | bcftools consensus -m /outdir/{accession}.mask.bed -H R -p {accession} /outdir/{accession}.vcf.gz > /outdir/{accession}.consensus.fa\''
            subprocess.check_call(consensus_cmd, shell=True)
    return is_plot_worthy

def run(accession_list, pe1, outdir, prefix, blast_db, read_length, pe2=None, plot_threshold=60, fa_threshold=70,threads=None):
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    blast_db_dir = os.path.abspath(os.path.dirname(blast_db))
    blast_db_name = os.path.basename(blast_db)

    pe1 = os.path.abspath(pe1)
    if pe2 is not None:
        pe2 = os.path.abspath(pe2)

    descriptions = {}
    descriptions_lock = threading.Lock()

    # Determine number of threads
    if threads is None or threads <= 0:
        threads = os.cpu_count()
        if threads is None:
            threads = 1

    print(f"Using a thread pool with {threads} workers.")

    # Create and run the thread pool
    plot_worthy_accessions = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_accession, acc, pe1, pe2, outdir, blast_db_dir, blast_db_name, read_length,
                                   plot_threshold, fa_threshold, docker, descriptions_lock, descriptions): acc for acc in accession_list}

        for future in futures:
            acc = futures[future]
            try:
                is_plot_worthy = future.result()
                if is_plot_worthy:
                    plot_worthy_accessions.append(acc)
            except Exception as e:
                print(f"Processing for accession {acc} failed: {e}")
    with open(f"{outdir}/{prefix}.final.txt", "w") as outfile1:
        num=0
        for accession in accession_list:
            cov_file = f"{outdir}/{accession}.cov.txt"
            if os.path.exists(cov_file):
                with open(cov_file, "r") as f_cov:
                    for line in f_cov:
                        line = line.strip()
                        if not line.startswith("#"):
                            array = line.split("\t")
                            reads = float(array[6]) + float(array[7])
                            if float(array[5]) > max(int(read_length) * 3, 500) or float(array[1]) >= 10 or float(array[9]) >= 10 or reads >=3:
                                outfile1.write(f"{descriptions.get(array[0], 'N/A')}\t{array[1]}\t{array[2]}\t{array[3]}\t{array[4]}\t{array[5]}\t{array[6]}\t{array[7]}\t{array[8]}\t{array[9]}\n")
                        else:
                            num += 1
                            if num==1:
                                outfile1.write(f"{line}\n")

    # Step 4: Plotting (single-threaded, after all data is collected)
    if plot_worthy_accessions:
        print("Step 4: Generating plots...")
        with open(f"{outdir}/{prefix}.depth.plot", "w") as depth_plot_file:
            for acc in plot_worthy_accessions:
                depth_cmd = f'docker run --rm -v {outdir}:/outdir {docker} sh -c \'export PATH=/opt/conda/bin:$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa /outdir/{acc}.bam >> /outdir/{prefix}.depth.plot\''
                subprocess.check_call(depth_cmd, shell=True)

        try:
            df = pd.read_csv(f"{outdir}/{prefix}.depth.plot", sep='\t', header=None, names=['Chr', 'Pos', 'Depth'])
        except pd.errors.EmptyDataError:
            print("The depth data file is empty. No plots will be generated.")
            return

        unique_chromosomes = df['Chr'].unique()

        # Part A: Plotting each chromosome individually
        print("Generating individual plots for each chromosome...")
        for chr_name in unique_chromosomes:
            plt.figure(figsize=(15, 6))
            chr_data = df[df['Chr'] == chr_name]
            plt.plot(chr_data['Pos'], chr_data['Depth'], color='blue', alpha=0.7)
            plt.yscale('log')
            plt.tick_params(axis='y', which='both', labelsize=12)
            median_depth = chr_data['Depth'].median()
            low_coverage_bases = (chr_data['Depth'] < 10).sum()
            plt.axhline(y=median_depth, color='r', linestyle='--', linewidth=2, label=f'Median: {int(median_depth)}')
            plt.axhline(y=10, color='gray', linestyle='--', linewidth=2, label=f'<10X: {low_coverage_bases} bp')
            title_desc = descriptions.get(chr_name, chr_name)
            plt.title(f"Chromosome: {chr_name}\nSpecies Description: {title_desc}", fontsize=16)
            plt.xlabel("Position (bp)", fontsize=12)
            plt.ylabel("Depth", fontsize=12)
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.legend(loc='lower right', title="", frameon=True)
            plt.tight_layout()
            plt.savefig(f"{outdir}/{prefix}_{chr_name}.png", dpi=300)
            plt.close()
            print(f"Individual plot saved for {chr_name}.")

        # Part B: Plotting all chromosomes as subplots
        print("\nGenerating a combined plot with subplots for all chromosomes...")
        num_plots = len(unique_chromosomes)
        if len(unique_chromosomes) <= 30:
            ncols = 2
            nrows = (num_plots + ncols - 1) // ncols

            labels = {chr_name: descriptions.get(chr_name, chr_name) for chr_name in unique_chromosomes}
            longest_label_len = max(len(label) for label in labels.values()) if labels else 0
            base_width = 15
            title_extra_width = max(0, (longest_label_len - 50) * 0.1)
            fig_width = ncols * (base_width + title_extra_width)
            fig_height = nrows * 6

            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width, fig_height), sharex=False)
            if num_plots == 1:
                axes = np.array([axes])
            axes = axes.flatten()

            for i, chr_name in enumerate(unique_chromosomes):
                ax = axes[i]
                chr_data = df[df['Chr'] == chr_name]
                ax.plot(chr_data['Pos'], chr_data['Depth'], color='blue', alpha=0.7)
                ax.set_yscale('log')
                median_depth = chr_data['Depth'].median()
                low_coverage_bases = (chr_data['Depth'] < 10).sum()
                ax.axhline(y=median_depth, color='r', linestyle='--', linewidth=2)
                ax.axhline(y=10, color='gray', linestyle='--', linewidth=2)
                legend_handles = [
                    plt.Line2D([0], [0], linestyle='None', marker='None', color='w'),
                    plt.Line2D([0], [0], linestyle='--', color='gray'),
                    plt.Line2D([0], [0], linestyle='--', color='r')
                ]
                legend_labels = [
                    "",  # An empty line for spacing
                    f'<10X: {low_coverage_bases} bp',
                    f'Median: {int(median_depth)}'
                ]
                ax.legend(legend_handles, legend_labels, loc='lower right', frameon=True)
                title = labels.get(chr_name, chr_name)
                ax.set_title(f"{title}", fontsize=14)
                ax.set_xlabel("Position (bp)", fontsize=10)
                ax.set_ylabel("Depth", fontsize=10)
                ax.grid(True, linestyle='--', alpha=0.6)

            for i in range(num_plots, len(axes)):
                fig.delaxes(axes[i])

            plt.tight_layout()
            plt.savefig(f"{outdir}/{prefix}_all_chromosomes_subplots.png", dpi=300)
            plt.close()
            print(f"Subplot plot for all chromosomes saved to {outdir}/{prefix}_all_chromosomes_subplots.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant calling + consensus with N mask and coverage plotting.")
    parser.add_argument("-p1", "--pe1", help="pe1 fastq", required=True)
    parser.add_argument("-p2", "--pe2", help="pe2 fastq", default=None)
    parser.add_argument("-p", "--prefix", help="output directory", required=True)
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-a", "--accession", required=True, nargs='+')
    parser.add_argument("-d", "--blast-db", help="Path to the local BLAST database for getting species descriptions",required=True)
    parser.add_argument("-plot", "--plot", default=60, type=int, help="plot when Covered_percent >default=60")
    parser.add_argument("-fa", "--fa", default=70, type=int,help="output consensus fasta when Covered_percent >default=70")
    parser.add_argument("-l", "--read_length", type=int, help="read length", required=True)
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count(), help="Number of threads to use (default: number of CPU cores)")
    args = parser.parse_args()
    run(args.accession, args.pe1, args.outdir, args.prefix, args.blast_db, args.read_length, args.pe2, args.plot,args.fa, args.threads)
