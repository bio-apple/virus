import subprocess
import os
import argparse

def run(bam,ref,outdir,prefix):
    bam=os.path.abspath(bam)
    out=os.path.abspath(outdir)+"/"+prefix
    result = subprocess.run(f"samtools view -H {bam}", shell=True, capture_output=True, text=True, check=True)
    if not "SO:coordinate" in result.stdout:
        sort_bam=os.path.splitext(bam)[0]+"sort.bam"
        subprocess.run(f"samtools sort -o {sort_bam} {bam}", shell=True, check=True)
    else:
        sort_bam=bam
    if not os.path.exists(sort_bam + ".bai"):
        subprocess.check_call(f'samtools index {sort_bam}',shell=True)

    cmd=(f'bcftools mpileup -Ou -f {ref} {sort_bam} | '
         f'bcftools call --ploidy 1 -mv -Oz -o {out}.calls.vcf.gz && '
         f'bcftools index {out}.calls.vcf.gz && '
         f'cat {ref} | bcftools consensus {out}.calls.vcf.gz -p {prefix} > {out}.consensus.fa')
    subprocess.check_call(cmd,shell=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser("Variant calling and calculation of the consensus sequence using bcftools.")
    parser.add_argument("-b","--bam",help="sort bam file",required=True)
    parser.add_argument("-r","--ref",help="reference fasta",required=True)
    parser.add_argument("-o","--outdir",help="output directory",required=True)
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    args=parser.parse_args()
    run(args.bam,args.ref,args.outdir,args.prefix)
