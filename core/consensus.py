import os
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics
import math

docker = "virus:latest"

def run(bam, ref, outdir, prefix):
    bam = os.path.abspath(bam)
    ref = os.path.abspath(ref)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    bam_in = f"/raw_data/{os.path.basename(bam)}"
    ref_in = f"/ref/{os.path.basename(ref)}"
    sort_bam = os.path.splitext(os.path.basename(bam))[0] + ".sort.bam"
    sort_bam_in = f"/outdir/{sort_bam}"
    depth_txt = f"/outdir/{prefix}.depth.txt"
    volume_args = f"-v {bam}:{bam_in} -v {ref}:{ref_in} -v {outdir}:/outdir"

    # Step 1: Check sort
    header_cmd = (
        f"docker run --rm {volume_args} {docker} "
        f"sh -c 'export PATH=/opt/conda/bin/:$PATH && samtools view -H {bam_in}'"
    )
    result = subprocess.run(header_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)

    if "SO:coordinate" not in result.stdout:
        sort_cmd = (
            f"docker run --rm {volume_args} {docker} "
            f"sh -c 'export PATH=/opt/conda/bin/:$PATH && samtools sort -o {sort_bam_in} {bam_in}'"
        )
        subprocess.check_call(sort_cmd, shell=True)
    else:
        sort_bam = os.path.basename(bam)
        sort_bam_in = bam_in

    # Step 2: Index
    bai_path = os.path.join(outdir, sort_bam + ".bai")
    if not os.path.exists(bai_path):
        index_cmd = (
            f"docker run --rm {volume_args} {docker} "
            f"sh -c 'export PATH=/opt/conda/bin/:$PATH && samtools index {sort_bam_in}'"
        )
        subprocess.check_call(index_cmd, shell=True)

    # Step 3: Coverage depth
    depth_cmd = (
        f"docker run --rm {volume_args} {docker} "
        f"sh -c 'export PATH=/opt/conda/bin/:$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa {sort_bam_in} > {depth_txt}'"
    )
    subprocess.check_call(depth_cmd, shell=True)

    # Step 4: Mask low coverage regions (<10x) and call consensus
    mask_bed = f"{outdir}/{prefix}.mask.bed"
    with open(f"{outdir}/{prefix}.depth.txt") as infile, open(mask_bed, "w") as outbed:
        for line in infile:
            chrom, pos, depth = line.strip().split("\t")
            if int(depth) < 10:
                outbed.write(f"{chrom}\t{int(pos)-1}\t{pos}\n")  # BED format

    consensus_cmd = (
        f"docker run --rm {volume_args} {docker} "
        f"sh -c 'export PATH=/opt/conda/bin/:$PATH && bcftools mpileup -Ou -f {ref_in} {sort_bam_in} | "
        f"bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.vcf.gz && "
        f"bcftools index /outdir/{prefix}.vcf.gz && "
        f"cat {ref_in} | bcftools consensus -m /outdir/{prefix}.mask.bed "
        f"-p {prefix} /outdir/{prefix}.vcf.gz > /outdir/{prefix}.consensus.fa'"
    )
    subprocess.check_call(consensus_cmd, shell=True)

    # Step 5: Plot coverage
    df = pd.read_csv(f"{outdir}/{prefix}.depth.txt", sep="\t", names=["ref", "pos", "depth"])
    chroms = df["ref"].unique()
    cols = 2
    rows = math.ceil(len(chroms) / cols)
    fig, axs = plt.subplots(rows, cols, figsize=(12, 4 * rows), squeeze=False)

    for idx, chrom in enumerate(chroms):
        sub_df = df[df["ref"] == chrom]
        median = statistics.median(sub_df["depth"])
        low_cov = (sub_df["depth"] < 10).sum()

        r, c = idx // cols, idx % cols
        ax = axs[r][c]
        ax.plot(sub_df["pos"], sub_df["depth"], linewidth=0.8)
        ax.axhline(10, linestyle="--", color="gray", label=f"<10X: {low_cov}bp")
        ax.axhline(median, linestyle="--", color="red", label=f"Median: {median:.0f}")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1, top=max(10, sub_df["depth"].max()) * 2)
        ax.set_title(f"Chromosome: {chrom}")
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Depth")
        ax.legend(fontsize=8)
        ax.margins(x=0.01)

    for idx in range(len(chroms), rows * cols):
        fig.delaxes(axs[idx // cols][idx % cols])

    plt.tight_layout()
    plt.savefig(f"{outdir}/{prefix}.coverage.png", dpi=300)
    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant calling + consensus with N mask and coverage plotting.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output")
    args = parser.parse_args()
    run(args.bam, args.ref, args.outdir, args.prefix)