import pandas as pd
import matplotlib.pyplot as plt
import subprocess, argparse
import os
import numpy as np

docker = "virus:latest"

def run(accession,pe1,outdir, prefix, blast_db,read_length,pe2=None,plot=60,fa=70):
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    blast_db_dir = os.path.abspath(os.path.dirname(blast_db))
    blast_db_name = blast_db_dir.split("/")[-1]
    pe1=os.path.abspath(pe1)
    if pe2 is not None:
        pe2=os.path.abspath(pe2)

    cmd = f'docker run --rm -v {outdir}:/outdir -v {blast_db_dir}:/ref {docker} sh -c \''

    outfile1=open(f"{outdir}/{prefix}.cov.txt","w")
    outfile2=open(f"{outdir}/{prefix}.final.txt","w")

    num,chr_plot = 0,[]
    for key in accession:
        num+=1
        fasta_cmd = cmd + (f"export PATH=/opt/conda/envs/kraken2/bin/:$PATH && blastdbcmd -db /ref/{blast_db_name} -entry {key} > /outdir/{key}.fasta\'")
        print(fasta_cmd)
        subprocess.check_call(fasta_cmd, shell=True)
        
        cov=""
        if pe2 is not None:
            cov = (f'docker run '
                   f'-v {pe1}:/raw_data/{pe1.split("/")[-1]} '
                   f'-v {pe2}:/raw_data/{pe2.split("/")[-1]} '
                   f'-v {outdir}:/outdir/ {docker} '
                   f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
                   f'bbwrap.sh ref=/outdir/{key}.fasta in1=/raw_data/{pe1.split("/")[-1]} in2=/raw_data/{pe2.split("/")[-1]} out=/outdir/{key}.sam.gz && '
                   f'pileup.sh in=/outdir/{key}.sam.gz out=/outdir/{key}.cov.txt\'')
        else:
            cov = (f'docker run '
                   f'-v {pe1}:/raw_data/{pe1.split("/")[-1]} '
                   f'-v {outdir}:/outdir/ {docker} '
                   f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
                   f'bbwrap.sh ref=/outdir/{key}.fasta in=/raw_data/{pe1.split("/")[-1]} out=/outdir/{key}.sam.gz && '
                   f'pileup.sh in=/outdir/{key}.sam.gz out=/outdir/{key}.cov.txt\'')
        print(cov)
        subprocess.check_call(cov, shell=True)
        with open(f"{outdir}/{key}.cov.txt","r") as f:
            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    array = line.split("\t")
                    if not (float(array[1]) == 0 and float(array[4]) == 0 and float(array[9]) == 0):
                        outfile1.write(f"{line}\n")#raw
                    if float(array[5]) > max(int(read_length) * 3, 500) or float(array[1]) >= 10 or float(array[9]) >= 10:
                        outfile2.write(f"{line}\n")#pos
                        if float(array[4]) >= plot:
                            chr_plot.append(array[0])
                            #refernce index
                            subprocess.check_call(
                                cmd+f'export PATH=/opt/conda/envs/kraken2/bin:/opt/conda/bin:$PATH && cd /outdir/ && '
                                f'bowtie2-build {key}.fasta {key}.fasta && '
                                f'samtools faidx {key}.fasta && '
                                f'bwa index -a bwtsw {key}.fasta\'', shell=True)

                            #mapping
                            mapping=f'docker run -v {pe1}:/raw_data/{pe1.split("/")[-1]} -v {outdir}:/outdir/ '
                            if pe2 is not None:
                                mapping+=f'-v {pe2}:/raw_data/{pe2.split("/")[-1]} '
                            mapping+=f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bowtie2 --threads 48 -x /outdir/{key}.fasta '
                            if pe2 is not None:
                                mapping += f'-1 /raw_data/{pe1.split("/")[-1]} -2 /raw_data/{pe2.split("/")[-1]}|samtools view -bh |samtools sort > /outdir/{key}.bam && samtools index /outdir/{key}.bam\''
                            else:
                                mapping+= f'-U {pe1}|samtools view -bh |samtools sort > /outdir/{key}.bam && samtools index /outdir/{key}.bam\''
                            print(mapping)
                            subprocess.check_call(mapping, shell=True)

                            #depth
                            depth_cmd = cmd + (f"export PATH=/opt/conda/bin:$PATH && "
                                               f"samtools depth -J -d 8000 -Q 0 -q 20 -aa /outdir/{key}.bam > /outdir/{key}.depth.txt\'")
                            subprocess.check_call(depth_cmd, shell=True)
                            #mask bed
                            with open(f"{outdir}/{key}.depth.txt", "r") as infile,open(f"{outdir}/{key}.mask.bed", "w") as outbed:
                                for line in infile:
                                    chrom, pos, depth = line.strip().split("\t")
                                    if int(depth) < 10:
                                        outbed.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")
                            subprocess.check_call(f'cat {outdir}/{key}.depth.txt >>{outdir}/{prefix}.depth.plot', shell=True)

                        if float(array[4]) >= fa:  # output fasta
                            print("variant calling...")
                            vcf_cmd = cmd + (f'export PATH=/opt/conda/bin/:$PATH && '
                                             f'bcftools mpileup -Ou -f /outdir/{key}.fasta /outdir/{key}.bam | '
                                             f'bcftools call --ploidy 1 -mv -Oz -o /outdir/{key}.vcf.gz\'')
                            print(vcf_cmd)
                            subprocess.check_call(vcf_cmd, shell=True)

                            print("consensus...")
                            consensus = cmd + (f'export PATH=/opt/conda/bin/:$PATH && '
                                               f'bcftools index /outdir/{key}.vcf.gz && cat /outdir/{key}.fasta | '
                                               f'bcftools consensus -m /outdir/{key}.mask.bed -H R -p {key} /outdir/{key}.vcf.gz > /outdir/{key}.consensus.fa\'')
                            print(consensus)
                            subprocess.check_call(consensus, shell=True)
                else:
                    if num==1:
                        outfile1.write(f"{line}\n")
                        outfile2.write(f"{line}\n")
        subprocess.check_call(f'cd {outdir} && rm -rf {key}.fasta* {key}.cov.txt {key}.sam.gz {key}.bam* {key}.depth.txt',shell=True)

    descriptions={}

    if len(chr_plot)!=0:
        print("Step 7: Generating plots...")

        # Load data for plotting
        try:
            df = pd.read_csv(f"{outdir}/{prefix}.depth.plot", sep='\t', header=None, names=['Chr', 'Pos', 'Depth'])
        except pd.errors.EmptyDataError:
            print("The depth data file is empty. No plots will be generated.")
            return

        # Get unique chromosomes from the DataFrame
        unique_chromosomes = df['Chr'].unique()


        # Part A: Plotting each chromosome individually and saving as separate files
        print("Generating individual plots for each chromosome...")
        for chr_name in unique_chromosomes:
            plt.figure(figsize=(15, 6))

            # Get data for the current chromosome
            chr_data = df[df['Chr'] == chr_name]

            # Plot the depth
            plt.plot(chr_data['Pos'], chr_data['Depth'], color='blue', alpha=0.7)

            # Set y-axis to a logarithmic scale
            plt.yscale('log')
            plt.tick_params(axis='y', which='both', labelsize=12)

            # Calculate key statistics for the legend
            median_depth = chr_data['Depth'].median()
            low_coverage_bases = (chr_data['Depth'] < 10).sum()

            # Add horizontal lines
            plt.axhline(y=median_depth, color='r', linestyle='--', linewidth=2, label=f'Median: {int(median_depth)}')
            plt.axhline(y=10, color='gray', linestyle='--', linewidth=2, label=f'<10X: {low_coverage_bases} bp')

            # Check if a description exists before adding it to the title
            title_desc = descriptions.get(chr_name)
            if title_desc:
                # If a description exists, include both chromosome and species info
                title_str = f"Chromosome: {chr_name}\nSpecies Description: {title_desc}"
            else:
                # If no description exists, only show the chromosome name
                title_str = f"Chromosome: {chr_name}"

            plt.title(title_str, fontsize=16)
            plt.xlabel("Position (bp)", fontsize=12)
            plt.ylabel("Depth", fontsize=12)
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.legend(loc='lower right', title="", frameon=True)

            plt.tight_layout()

            # Save the individual plot
            plt.savefig(f"{outdir}/{prefix}_{chr_name}.png", dpi=300)
            plt.close()
            print(f"Individual plot saved for {chr_name}.")

        # Part B: Plotting all chromosomes as subplots on a single large figure
        print("\nGenerating a combined plot with subplots for all chromosomes...")
        num_plots = len(unique_chromosomes)
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

            # Calculate key statistics for the legend
            median_depth = chr_data['Depth'].median()
            low_coverage_bases = (chr_data['Depth'] < 10).sum()

            # Add horizontal lines without labels
            ax.axhline(y=median_depth, color='r', linestyle='--', linewidth=2)
            ax.axhline(y=10, color='gray', linestyle='--', linewidth=2)

            # Create custom legend for the statistics
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
    parser.add_argument("-p1","--pe1", help="pe1 fastq",required=True)
    parser.add_argument("-p2","--pe2",help="pe2 fastq",default=None)
    parser.add_argument("-p","--prefix", help="output directory",required=True)
    parser.add_argument("-o","--outdir", required=True, help="Output directory")
    parser.add_argument("-a", "--accession", required=True, nargs='+')
    parser.add_argument("-d","--blast-db", default=None,help="Path to the local BLAST database for getting species descriptions (optional).")
    parser.add_argument("-plot","--plot", default=60, type=int, help="plot when Covered_percent >default=60")
    parser.add_argument("-fa","--fa", default=70, type=int, help="output consensus fasta when Covered_percent >default=70")
    parser.add_argument("-l", "--read_length", type=int, help="read length",required=True)
    args = parser.parse_args()
    run(args.accession, args.pe1, args.outdir,args.prefix, args.blast_db,args.read_length,args.plot, args.fa,args.pe2)