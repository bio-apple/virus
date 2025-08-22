import pandas as pd
import matplotlib.pyplot as plt
import subprocess, argparse
import os
import numpy as np
import random

docker = "virus:latest"

def run(bam, outdir, prefix, blast_db,read_length):
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    bam = os.path.abspath(bam)
    bam_file_name=bam.split("/")[-1]
    blast_db_dir = os.path.abspath(os.path.dirname(blast_db))
    blast_db_name = blast_db_dir.split("/")[-1]
    cmd = f'docker run --rm -v {bam}:/raw_data/{bam_file_name} -v {outdir}:/outdir -v {blast_db_dir}:/ref {docker} sh -c \''

    print("Step 1: Calculating sequencing depth and coverage summary...")
    depth_cmd =cmd+ (f"export PATH=/opt/conda/bin:$PATH && "
           f"samtools depth -J -d 8000 -Q 0 -q 20 -aa /raw_data/{os.path.basename(bam)} > /outdir/{prefix}.depth.txt && "
           f"pileup.sh in=/raw_data/{os.path.basename(bam)} out=/outdir/{prefix}.cov\'")
    subprocess.check_call(depth_cmd, shell=True)

    print("Step 2: Filtering chromosomes based on coverage...")
    chr_raw,chr_pos,chr_plot = [],[],[]
    positive,plot=10,70
    with open(f"{outdir}/{prefix}.cov", "r") as infile:
        for line in infile:
            line = line.strip()
            if not line.startswith("#"):
                array = line.split("\t")
                if not (float(array[1]) == 0 and float(array[4]) == 0 and float(array[9]) == 0):
                    chr_raw.append(array[0])
                if float(array[5]) >max(int(read_length)*3,500) or float(array[1])>=10 or float(array[9])>=10:
                    chr_pos.append(array[0])
                    if float(array[4]) >= plot:#plot
                        chr_plot.append(array[0])

    print("Step 3: Extracting depth data for filtered chromosomes...")
    mask_bed = f"{outdir}/{prefix}.mask.bed"
    with open(f"{outdir}/{prefix}.depth.txt", "r") as infile,open(f"{outdir}/{prefix}.depth.plot", "w") as outfile,open(mask_bed, "w") as outbed:
        for line in infile:
            chrom, pos, depth = line.strip().split("\t")
            if chrom in chr_plot:
                outfile.write(line)
                if int(depth) < positive:
                    outbed.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")

    descriptions = {}
    if blast_db:
        print("Step 4: Fetching descriptions and fasta sequence from BLAST database...")
        if os.path.exists(f"{outdir}/ref.final.fasta"):
            os.remove(f"{outdir}/ref.final.fasta")
        for key in chr_raw:#raw
            blast_cmd =cmd+f"export PATH=/opt/conda/envs/kraken2/bin/:$PATH && blastdbcmd -db /ref/{blast_db_name} -entry {key} -outfmt \"%a|%t\"\'"
            process = subprocess.Popen(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stdout, stderr = process.communicate()
            if process.returncode == 0 and stdout.strip():
                all_descriptions = [line.strip() for line in stdout.strip().split('\n') if line.strip()]
                if all_descriptions:
                    chosen_description = random.choice(all_descriptions)
                    descriptions[key] = chosen_description
            if key in chr_plot:#plot
                fasta_cmd = cmd+(f"export PATH=/opt/conda/envs/kraken2/bin/:$PATH && "
                             f"blastdbcmd -db /ref/{blast_db_name} -entry {key} >> /outdir/ref.final.fasta\'")
                subprocess.check_call(fasta_cmd, shell=True)

        print("Step 5: variant calling...")
        vcf_cmd = cmd+(f'export PATH=/opt/conda/bin/:$PATH && '
                        f'bcftools mpileup -Ou -f /outdir/ref.final.fasta /raw_data/{bam_file_name} | '
                        f'bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.vcf.gz\'')
        print(vcf_cmd)
        subprocess.check_call(vcf_cmd, shell=True)

        print("Step 6: consensus...")
        consensus = cmd + (f'export PATH=/opt/conda/bin/:$PATH && '
                      f'bcftools index /outdir/{prefix}.vcf.gz && cat /outdir/ref.final.fasta | '
                      f'bcftools consensus -m /outdir/{prefix}.mask.bed -p {prefix} /outdir/{prefix}.vcf.gz > /outdir/{prefix}.consensus.fa\'')
        print(consensus)
        subprocess.check_call(consensus, shell=True)

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
    with open(f"{outdir}/{prefix}.cov", "r") as infile, \
            open(f"{outdir}/{prefix}.cov.txt", "w") as outfile, \
            open(f"{outdir}/{prefix}.final.cov.txt", "w") as outfile1:
        for line in infile:
            line = line.strip()
            if line.startswith("#"):
                outfile.write(line + "\n")
                outfile1.write(line + "\n")
                continue

            array = line.split("\t")

            # Write to outfile if the chromosome is in chr_list_raw
            if array[0] in chr_raw:
                if array[0] in descriptions:
                    outfile.write(f"{descriptions[array[0]]}")
                    for i in range(1, len(array)):
                        outfile.write(f"\t{array[i]}")
                    outfile.write("\n")
                else:
                    outfile.write(line + "\n")

            # Write to outfile1 if the chromosome is in chr_list
            if array[0] in chr_pos:
                if array[0] in descriptions:
                    outfile1.write(f"{descriptions[array[0]]}")
                    for i in range(1, len(array)):
                        outfile1.write(f"\t{array[i]}")
                    outfile1.write("\n")
                else:
                    outfile1.write(line + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant calling + consensus with N mask and coverage plotting.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output")
    parser.add_argument("-d", "--blast-db", default=None,help="Path to the local BLAST database for getting species descriptions (optional).")
    parser.add_argument("-l", "--read_length", type=int, help="read length",required=True)
    args = parser.parse_args()
    run(args.bam, args.outdir, args.prefix, args.blast_db,args.read_length)