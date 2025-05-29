import os
import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import statistics
import math

docker = "virus:latest"

def run(bam, outdir, prefix, ref=None,chrom_list=None):
    bam = os.path.abspath(bam)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    cmd=f'docker run --rm -v {bam}:/raw_data/{bam.split("/")[-1]} -v {outdir}:/outdir '
    # Step 1: Generate depth file
    depth_cmd = cmd + f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && samtools depth -J -d 8000 -Q 0 -q 20 -aa /raw_data/{bam.split("/")[-1]} > /outdir/{prefix}.depth.txt\''
    subprocess.check_call(depth_cmd, shell=True)
    if ref:
        ref=os.path.abspath(ref)
        cmd+=f'-v {ref}:/ref/{ref.split("/")[-1]} '
        # Step 2: Create mask bed for positions with depth < 10
        mask_bed = f"{outdir}/{prefix}.mask.bed"
        with open(f"{outdir}/{prefix}.depth.txt") as infile, open(mask_bed, "w") as outbed:
            for line in infile:
                chrom, pos, depth = line.strip().split("\t")
                if int(depth) < 10:
                    outbed.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")  # BED format zero-based

        # Step 3: Call variants and generate consensus with masked low coverage regions
        consensus_cmd = cmd+f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bcftools mpileup -Ou -f /ref/{ref.split("/")[-1]} /raw_data/{bam.split("/")[-1]} | bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.vcf.gz\''
        print(consensus_cmd)
        subprocess.check_call(consensus_cmd, shell=True)

        temp=cmd+f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bcftools index /outdir/{prefix}.vcf.gz && cat /ref/{ref.split("/")[-1]} | bcftools consensus -m /outdir/{prefix}.mask.bed -p {prefix} /outdir/{prefix}.vcf.gz > /outdir/{prefix}.consensus.fa\''
        print(temp)
        subprocess.check_call(temp, shell=True)

    # Step 4: Plot coverage
    df = pd.read_csv(f"{outdir}/{prefix}.depth.txt", sep="\t", names=["ref", "pos", "depth"])
    available_chroms = df["ref"].unique()
    chroms = available_chroms
    print(chroms)
    if chrom_list:
        chroms = [c for c in chrom_list.split(" ") if c in available_chroms]
        print("Selected chromosomes for plotting:", chroms)

    cols = 2
    rows = math.ceil(len(chroms) / cols)

    # 单个图的特殊处理：返回不是数组而是一个单个 ax
    if len(chroms) == 1:
        fig, ax = plt.subplots(figsize=(8, 4))
        axs = [[ax]]
    else:
        fig, axs = plt.subplots(rows, cols, figsize=(12, 4 * rows), squeeze=False)

    for idx, chrom in enumerate(chroms):
        sub_df = df[df["ref"] == chrom]
        median = statistics.median(sub_df["depth"])
        low_cov = (sub_df["depth"] < 10).sum()

        r, c = idx // cols, idx % cols
        ax = axs[r][c]
        ax.plot(sub_df["pos"], sub_df["depth"], linewidth=0.8)
        ax.axhline(10, linestyle="--", color="gray", label=f"<10X: {low_cov} bp")
        ax.axhline(median, linestyle="--", color="red", label=f"Median: {median:.0f}")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1, top=max(10, sub_df["depth"].max()) * 2)
        ax.set_title(f"Chromosome: {chrom}")
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Depth")
        ax.legend(fontsize=8)
        ax.margins(x=0.01)

    # 删除空的 subplot（仅在多个图时需要）
    if len(chroms) > 1:
        for idx in range(len(chroms), rows * cols):
            fig.delaxes(axs[idx // cols][idx % cols])

    plt.tight_layout()
    plt.savefig(f"{outdir}/{prefix}.coverage.png", dpi=300)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Variant calling + consensus with N mask and coverage plotting.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-r", "--ref", default=None, help="Reference FASTA (optional)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for output")
    parser.add_argument("-c", "--chrom-list",nargs='+',default=None, help="List of chromosomes to plot (default: all)")
    args = parser.parse_args()
    run(args.bam, args.outdir,args.prefix, args.ref, args.chrom_list)
