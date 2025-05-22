import os,sys,re
import subprocess
import argparse

# align reads with bowtie2 and sort bam with samtools
def run(ref_index_dir,outdir,prefix,R1,R2=None):
    outdir=os.path.abspath(outdir)
    if not os.path.exists(outdir):
        subprocess.check_call(f'mkdir -p {outdir}',shell=True)
    out=outdir+"/"+prefix
    ref_index=""
    for i in os.listdir(ref_index_dir):
        if i.endswith(".rev.2.bt2"):
            ref_index = i.split(".rev.2.bt2")[0]
    if R2 != None:
        cmd= f"bowtie2 --threads 48 -x {ref_index_dir}/{ref_index} -1 {R1} -2 {R2}|samtools view -bh |samtools sort > {out}.bam && samtools index /outdir/{out}.bam"
    else:
        cmd= f"bowtie2 --threads 48 -x {ref_index_dir}/{ref_index} -U {R1}|samtools view -bh |samtools sort > {out}.bam && samtools index {out}.bam"
    subprocess.check_call(cmd, shell=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser("Mapping reference.")
    parser.add_argument("-1","--R1",help="R1 fastq",required=True)
    parser.add_argument("-2","--R2",help="R2 fastq",required=True)
    parser.add_argument("-r","--ref",help="directory of bowtie2 reference index",required=True)
    parser.add_argument("-o","--outdir",help="output directory",default=os.getcwd())
    parser.add_argument("-p","--prefix",help="prefix of output")
    args=parser.parse_args()
    run(args.ref,args.outdir,args.prefix,args.R1,args.R2)
