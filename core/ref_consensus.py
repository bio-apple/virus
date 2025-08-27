import os
import argparse
import subprocess, statistics
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

docker = "virus:latest"


def run(bam, outdir, prefix, ref):
    bam = os.path.abspath(bam)
    bam_name = bam.split("/")[-1]
    ref = os.path.abspath(ref)
    ref_name = ref.split("/")[-1]
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    cmd = f'docker run --rm -v {bam}:/raw_data/{bam_name} -v {outdir}:/outdir -v {os.path.dirname(ref)}:/ref/ {docker} sh -c \''
    depth_cmd = cmd + f'export PATH=/opt/conda/bin:$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa /raw_data/{bam_name} > /outdir/{prefix}.depth.txt\''
    subprocess.check_call(depth_cmd, shell=True)
    mask_bed = f"{outdir}/{prefix}.mask.bed"
    with open(f"{outdir}/{prefix}.depth.txt") as infile, open(mask_bed, "w") as outbed:
        for line in infile:
            chrom, pos, depth = line.strip().split("\t")
            if int(depth) < 10:
                outbed.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")

    vcf_cmd = cmd + (f'export PATH=/opt/conda/bin/:$PATH && '
                     f'bcftools mpileup -Ou -f /ref/{ref_name} /raw_data/{bam_name} | '
                     f'bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.vcf.gz\'')
    print(vcf_cmd)
    subprocess.check_call(vcf_cmd, shell=True)

    consensus_fa_path = f"{outdir}/{prefix}.consensus.fa"
    temp = cmd + (f'export PATH=/opt/conda/bin/:$PATH && '
                  f'bcftools index /outdir/{prefix}.vcf.gz && cat /ref/{ref_name} | '
                  f'bcftools consensus -m /outdir/{prefix}.mask.bed -H R -p {prefix} /outdir/{prefix}.vcf.gz > /outdir/{os.path.basename(consensus_fa_path)}\'')
    print(temp)
    subprocess.check_call(temp, shell=True)

    ##############plot coverage###########################
    df = pd.read_csv(f"{outdir}/{prefix}.depth.txt",
                     sep='\t',
                     names=["ref", "pos", "depth"]
                     )

    # Get unique chromosomes
    chromosomes = df['ref'].unique()

    # Loop through each chromosome and generate a separate plot
    for chrom in chromosomes:
        df_chrom = df[df['ref'] == chrom].copy()
        median_depth = statistics.median(df_chrom["depth"])
        low_coverage_bases = (df_chrom['depth'] < 10).sum()
        plt.figure(figsize=[10, 4])
        plt.axhline(median_depth, linestyle='--', color='red', linewidth=1, label="median: %.0f" % median_depth)
        plt.axhline(10, linestyle='--', color='grey', linewidth=1, label="<10X(%s bp)" % (low_coverage_bases))

        max_val = np.max(df_chrom["depth"])
        maxlog10 = np.ceil(np.log10(max_val)) if max_val > 1 else 1
        plt.ylim(top=10 ** maxlog10)

        plt.title(f"Sample: {prefix}\nChromosome: {chrom}", fontsize=10, wrap=True)
        plt.xlabel("Position along genome [bp]")
        plt.ylabel("Coverage depth")
        plt.yscale("log")
        plt.margins(x=0.01)
        plt.legend()
        plt.ylim(bottom=1)
        plt.plot(df_chrom["pos"], df_chrom["depth"])
        plt.savefig(f"{outdir}/{prefix}.{chrom}.coverage.png", dpi=300)
        plt.close()  # Close the figure to free up memory


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant calling + consensus with N mask and coverage plotting.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output")
    parser.add_argument("-r", "--ref", required=True, help="Reference fasta file")
    args = parser.parse_args()
    run(args.bam, args.outdir, args.prefix, args.ref)