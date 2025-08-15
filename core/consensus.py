import os
import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import statistics
import math
import random

docker = "virus:latest"
# Unified command prefix to set the PATH inside the Docker container
# This ensures all tools are found in the specified Conda environment
DOCKER_CMD_PREFIX = "export PATH=/opt/conda/envs/kraken2/bin/:$PATH && "


def get_chrom_descriptions(chroms, blast_db):
    """
    Constructs and runs a docker command with blastdbcmd to get descriptions.
    If multiple descriptions are found for a chromosome, one is chosen at random.
    """
    descriptions = {}
    if not blast_db:
        print("No local BLAST database provided. Skipping species description retrieval.")
        return descriptions

    print("Using blastdbcmd to extract chromosome descriptions from the local database...")

    # We need to map the local blast_db path to a path inside the Docker container
    db_dir = os.path.dirname(os.path.abspath(blast_db))
    db_basename = os.path.basename(blast_db)

    # Docker mount point for the blastdb directory
    docker_db_dir = "/blast_db"

    # Build the docker command to run blastdbcmd for each chromosome
    for chrom in chroms:
        # Construct and run the blastdbcmd command inside docker
        # The DOCKER_CMD_PREFIX ensures blastdbcmd is found
        cmd = f'docker run --rm -v {db_dir}:{docker_db_dir} {docker} sh -c "{DOCKER_CMD_PREFIX} blastdbcmd -db {docker_db_dir}/{db_basename} -entry {chrom} -outfmt %t"'

        try:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       universal_newlines=True)
            stdout, stderr = process.communicate()

            if process.returncode == 0 and stdout.strip():
                # Split the output by newline to get all descriptions
                all_descriptions = [line.strip() for line in stdout.strip().split('\n') if line.strip()]

                if all_descriptions:
                    # Randomly select one description
                    chosen_description = random.choice(all_descriptions)
                    descriptions[chrom] = chosen_description
                else:
                    print(f"Warning: Could not find a description for chromosome '{chrom}' in the database.")
            elif stderr:
                pass
            else:
                print(f"Warning: Could not find a description for chromosome '{chrom}' in the database.")
        except FileNotFoundError:
            print("Error: The blastdbcmd command was not found. Please ensure BLAST+ is installed and in your PATH.")
            return {}
        except Exception as e:
            print(f"An error occurred while processing blastdbcmd output: {e}")
            return {}

    return descriptions


def plot_coverage(df, chroms, outdir, prefix, suffix, descriptions):
    """
    Helper function to plot coverage for a given list of chromosomes.
    If there are 5 or fewer chromosomes, each is plotted in a separate file.
    Otherwise, a single multi-panel plot is generated.
    """
    if not chroms:
        print(f"No {suffix} chromosomes to plot. Skipping.")
        return

    # --- Plotting individual chromosomes if the count is low ---
    if len(chroms) <= 5:
        print(f"Found {len(chroms)} {suffix} chromosomes. Generating individual plots.")
        for chrom in chroms:
            sub_df = df[df["ref"] == chrom]
            if sub_df.empty:
                continue

            median = statistics.median(sub_df["depth"])
            low_cov = (sub_df["depth"] < 10).sum()

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(sub_df["pos"], sub_df["depth"], linewidth=0.8)
            ax.axhline(10, linestyle="--", color="gray", label=f"<10X: {low_cov} bp")
            ax.axhline(median, linestyle="--", color="red", label=f"Median: {median:.0f}")
            ax.set_yscale("log")
            ax.set_ylim(bottom=1, top=max(10, sub_df["depth"].max()) * 2)

            title_text = f"Chromosome: {chrom}"
            if descriptions and chrom in descriptions:
                desc = descriptions[chrom]
                title_text += f"\nSpecies Description: {desc}"

            ax.set_title(title_text, fontsize=12, wrap=True)
            ax.set_xlabel("Position (bp)")
            ax.set_ylabel("Depth")
            ax.legend()
            ax.margins(x=0.01)

            # Save the individual plot with the chromosome name in the filename
            plt.tight_layout()
            plt.savefig(f"{outdir}/{prefix}.{suffix}.{chrom}.png", dpi=300)
            plt.close(fig)
        return

    # --- Original logic for multi-panel plotting when count is high ---
    print(f"Found {len(chroms)} {suffix} chromosomes. Generating a single multi-panel plot.")

    # User's request was for 4, but the original code was 2. Setting it to 4
    cols = 4
    rows = math.ceil(len(chroms) / cols)

    # Base figure width for a 2-column layout
    base_fig_width = 12
    # Base height for a subplot with a single-line title
    base_subplot_height = 5.0

    # --- Dynamic Width Calculation ---
    max_desc_len = 0
    if descriptions:
        for desc in descriptions.values():
            # Consider both single lines and multi-line descriptions
            lines = desc.split('\n')
            for line in lines:
                max_desc_len = max(max_desc_len, len(line))

    # A heuristic to estimate required width.
    char_per_inch = 15  # A simple guess of characters per inch
    required_width_per_plot = max(base_fig_width / cols, max_desc_len / char_per_inch)

    # The final figure width is the required width per plot multiplied by the number of columns
    dynamic_fig_width = required_width_per_plot * cols

    # --- Dynamic Height Calculation (retained) ---
    row_heights = []
    for r in range(rows):
        start_idx = r * cols
        end_idx = min(start_idx + cols, len(chroms))

        max_height_in_row = base_subplot_height
        for idx in range(start_idx, end_idx):
            chrom = chroms[idx]
            if descriptions and chrom in descriptions:
                desc = descriptions[chrom]
                # Re-estimate lines based on the new dynamic width
                estimated_lines = desc.count('\n') + 1 + math.ceil(
                    len(desc) / (char_per_inch * required_width_per_plot))
                max_height_in_row = max(max_height_in_row, base_subplot_height + (estimated_lines * 0.5))

        row_heights.append(max_height_in_row)

    total_fig_height = sum(row_heights)

    # Create the figure with the dynamically calculated width and height
    fig, axs = plt.subplots(rows, cols, figsize=(dynamic_fig_width, total_fig_height), squeeze=False)

    for idx, chrom in enumerate(chroms):
        sub_df = df[df["ref"] == chrom]
        # Skip if no data for this chromosome
        if sub_df.empty:
            continue

        median = statistics.median(sub_df["depth"])
        low_cov = (sub_df["depth"] < 10).sum()

        r, c = idx // cols, idx % cols
        ax = axs[r][c]
        ax.plot(sub_df["pos"], sub_df["depth"], linewidth=0.8)
        ax.axhline(10, linestyle="--", color="gray", label=f"<10X: {low_cov} bp")
        ax.axhline(median, linestyle="--", color="red", label=f"Median: {median:.0f}")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1, top=max(10, sub_df["depth"].max()) * 2)

        # Update title with species description if available
        title_text = f"Chromosome: {chrom}"
        if descriptions and chrom in descriptions:
            desc = descriptions[chrom]
            title_text += f"\nSpecies Description: {desc}"

        ax.set_title(title_text, fontsize=10, loc='center', wrap=True)

        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Depth")
        ax.legend(fontsize=8)
        ax.margins(x=0.01)

    if len(chroms) > 1:
        for idx in range(len(chroms), rows * cols):
            fig.delaxes(axs[idx // cols][idx % cols])

    plt.suptitle(f"Coverage Plot - {suffix.replace('_', ' ').title()}", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f"{outdir}/{prefix}.{suffix}.png", dpi=300)
    plt.close(fig)


def run(bam, outdir, prefix, ref=None, chrom_list=None, blast_db=None):
    bam = os.path.abspath(bam)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Step 1: Extract chromosome names from BAM header
    print("Extracting chromosome names from BAM file...")
    try:
        bam_path = f"/raw_data/{os.path.basename(bam)}"
        cmd_get_chroms = f'docker run --rm -v {bam}:/raw_data/{os.path.basename(bam)} {docker} sh -c "export PATH=/opt/conda/bin:$PATH && samtools view -H {bam_path}"'
        header_output = subprocess.check_output(cmd_get_chroms, shell=True, universal_newlines=True)
        # Use regex to find all "SN:chromosome_name" lines
        all_chroms_from_bam = [line.split('\t')[1][3:] for line in header_output.split('\n') if line.startswith('@SQ')]
        if not all_chroms_from_bam:
            raise ValueError(
                "Could not extract any chromosome names from the BAM file. Please check if the file is valid.")
        print(f"Extracted {len(all_chroms_from_bam)} chromosomes.")
    except Exception as e:
        print(f"An error occurred while extracting chromosome names from the BAM file: {e}")
        return

    # Determine which chromosomes to plot based on user input
    chroms_to_process = all_chroms_from_bam
    if chrom_list:
        chroms_to_process = [c for c in chrom_list if c in all_chroms_from_bam]
        print(f"User specified {len(chroms_to_process)} chromosomes to process.")
    else:
        print("No chromosomes specified by the user; processing all chromosomes.")

    cmd = f'docker run --rm -v {bam}:/raw_data/{os.path.basename(bam)} -v {outdir}:/outdir '

    # Step 2: Generate depth file
    # The DOCKER_CMD_PREFIX is now part of the command string
    depth_cmd = cmd + f'{docker} sh -c \'export PATH=/opt/conda/bin:$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa /raw_data/{os.path.basename(bam)} > /outdir/{prefix}.depth.txt\''
    subprocess.check_call(depth_cmd, shell=True)

    if ref:
        ref = os.path.abspath(ref)
        cmd += f'-v {ref}:/ref/{os.path.basename(ref)} '
        # Step 3: Create mask bed for positions with depth < 10
        mask_bed = f"{outdir}/{prefix}.mask.bed"
        with open(f"{outdir}/{prefix}.depth.txt") as infile, open(mask_bed, "w") as outbed:
            for line in infile:
                chrom, pos, depth = line.strip().split("\t")
                if int(depth) < 10:
                    outbed.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")  # BED format zero-based

        # Step 4: Call variants and generate consensus with masked low coverage regions
        # The DOCKER_CMD_PREFIX is now part of the command string
        consensus_cmd = cmd + f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bcftools mpileup -Ou -f /ref/{os.path.basename(ref)} /raw_data/{os.path.basename(bam)} | export PATH=/opt/conda/bin/:$PATH && bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.vcf.gz\''
        print(consensus_cmd)
        subprocess.check_call(consensus_cmd, shell=True)

        temp = cmd + f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bcftools index /outdir/{prefix}.vcf.gz && cat /ref/{os.path.basename(ref)} | export PATH=/opt/conda/bin/:$PATH && bcftools consensus -m /outdir/{prefix}.mask.bed -p {prefix} /outdir/{prefix}.vcf.gz > /outdir/{prefix}.consensus.fa\''
        print(temp)
        subprocess.check_call(temp, shell=True)

    # Step 5: Plot coverage
    df = pd.read_csv(f"{outdir}/{prefix}.depth.txt", sep="\t", names=["ref", "pos", "depth"])

    high_coverage_chroms_smooth = []
    high_coverage_chroms_volatile = []
    low_coverage_chroms = []

    # Define the volatility threshold. This value may need to be adjusted based on your data.
    VOLATILITY_THRESHOLD = 0.5

    # Define the percentage of sequence to trim from each end (e.g., 1%)
    TRIM_PERCENT = 0.01

    print("Classifying chromosomes based on coverage volatility, ignoring ends...")

    for chrom in chroms_to_process:
        sub_df = df[df["ref"] == chrom]
        if sub_df.empty:
            continue

        median = statistics.median(sub_df["depth"])

        if median >= 10:
            # Get the total length of the chromosome
            chrom_length = sub_df["pos"].max()
            if chrom_length is None or chrom_length < 100:
                # If length is too short to trim, handle gracefully
                high_coverage_chroms_smooth.append(chrom)
                continue

            # Calculate the trim start and end positions
            trim_start = chrom_length * TRIM_PERCENT
            trim_end = chrom_length * (1 - TRIM_PERCENT)

            # Slice the DataFrame to exclude the ends
            trimmed_df = sub_df[(sub_df["pos"] > trim_start) & (sub_df["pos"] < trim_end)]

            if trimmed_df.empty:
                # If trimming results in no data, handle gracefully
                high_coverage_chroms_smooth.append(chrom)
                continue

            # Calculate volatility metric on the trimmed data: standard deviation of log-scaled depths
            # Filter out depths of 0 to avoid log10(0)
            log_depths = trimmed_df[trimmed_df["depth"] > 0]["depth"].apply(lambda x: math.log10(x))

            # Check if there are enough valid data points to compute standard deviation
            if len(log_depths) > 1:
                volatility = statistics.stdev(log_depths)

                if volatility < VOLATILITY_THRESHOLD:
                    high_coverage_chroms_smooth.append(chrom)
                else:
                    high_coverage_chroms_volatile.append(chrom)
            else:
                # If too few data points to calculate stdev, classify as smooth
                high_coverage_chroms_smooth.append(chrom)
        else:
            low_coverage_chroms.append(chrom)

    # Get descriptions for each group
    descriptions_all = get_chrom_descriptions(chroms_to_process, blast_db)

    descriptions_smooth = {chrom: descriptions_all.get(chrom, "") for chrom in high_coverage_chroms_smooth}
    descriptions_volatile = {chrom: descriptions_all.get(chrom, "") for chrom in high_coverage_chroms_volatile}
    descriptions_low_cov = {chrom: descriptions_all.get(chrom, "") for chrom in low_coverage_chroms}

    # Plot high coverage smooth chromosomes
    print("Plotting high-coverage 'smooth' chromosomes...")
    plot_coverage(df, high_coverage_chroms_smooth, outdir, prefix, "high_coverage_smooth", descriptions_smooth)

    # Plot high coverage volatile chromosomes
    print("Plotting high-coverage 'volatile' chromosomes...")
    plot_coverage(df, high_coverage_chroms_volatile, outdir, prefix, "high_coverage_volatile", descriptions_volatile)

    # Plot low coverage chromosomes
    print("Plotting low-coverage chromosomes...")
    plot_coverage(df, low_coverage_chroms, outdir, prefix, "low_coverage", descriptions_low_cov)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant calling + consensus with N mask and coverage plotting.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-r", "--ref", default=None, help="Reference FASTA (optional)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output")
    parser.add_argument("-c", "--chrom-list", nargs='+', default=None,
                        help="List of chromosomes to plot (default: all)")
    parser.add_argument("-d", "--blast-db", default=None,
                        help="Path to the local BLAST database for getting species descriptions (optional).")
    args = parser.parse_args()
    run(args.bam, args.outdir, args.prefix, args.ref, args.chrom_list, args.blast_db)