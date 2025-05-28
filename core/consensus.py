import subprocess
import os
import argparse

docker="fanyucai1/virus:latest"

def run(bam, ref, outdir, prefix):
    bam = os.path.abspath(bam)
    ref = os.path.abspath(ref)
    outdir = os.path.abspath(outdir)

    bam_in_container = f"/raw_data/{os.path.basename(bam)}"
    ref_in_container = f"/ref/{os.path.basename(ref)}"
    sort_bam_name = os.path.splitext(os.path.basename(bam))[0] + ".sort.bam"
    sort_bam_in_container = f"/outdir/{sort_bam_name}"

    # Docker volume mounts
    volume_args = (
        f"-v {bam}:{bam_in_container} "
        f"-v {ref}:{ref_in_container} "
        f"-v {outdir}:/outdir"
    )

    # 检查是否排序
    header_cmd = (
        f"docker run --rm {volume_args} {docker} "
        f"sh -c 'export PATH=/opt/conda/bin:$PATH && samtools view -H {bam_in_container}'"
    )
    result = subprocess.run(header_cmd, shell=True, capture_output=True, text=True, check=True)

    # 排序 BAM（如果未 coordinate）
    if "SO:coordinate" not in result.stdout:
        sort_cmd = (
            f"docker run --rm {volume_args} {docker} "
            f"sh -c 'export PATH=/opt/conda/bin:$PATH && "
            f"samtools sort -o {sort_bam_in_container} {bam_in_container}'"
        )
        subprocess.check_call(sort_cmd, shell=True)
    else:
        sort_bam_name = os.path.basename(bam)
        sort_bam_in_container = bam_in_container

    # 创建 BAM 索引
    bai_path = os.path.join(outdir, sort_bam_name + ".bai")
    if not os.path.exists(bai_path):
        index_cmd = (
            f"docker run --rm {volume_args} {docker} "
            f"sh -c 'export PATH=/opt/conda/bin:$PATH && "
            f"samtools index {sort_bam_in_container}'"
        )
        subprocess.check_call(index_cmd, shell=True)

    # 调用 bcftools
    bcftools_cmd = (
        f"docker run --rm {volume_args} {docker} "
        f"sh -c 'export PATH=/opt/conda/bin:$PATH && "
        f"bcftools mpileup -Ou -f {ref_in_container} {sort_bam_in_container} | "
        f"bcftools call --ploidy 1 -mv -Oz -o /outdir/{prefix}.calls.vcf.gz && "
        f"bcftools index /outdir/{prefix}.calls.vcf.gz && "
        f"cat {ref_in_container} | "
        f"bcftools consensus /outdir/{prefix}.calls.vcf.gz -p {prefix} > /outdir/{prefix}.consensus.fa'"
    )
    subprocess.check_call(bcftools_cmd, shell=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser("Variant calling and calculation of the consensus sequence using bcftools.")
    parser.add_argument("-b","--bam",help="sort bam file",required=True)
    parser.add_argument("-r","--ref",help="reference fasta",required=True)
    parser.add_argument("-o","--outdir",help="output directory",required=True)
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    args=parser.parse_args()
    run(args.bam,args.ref,args.outdir,args.prefix)
