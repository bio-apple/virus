import os,sys,re
import subprocess
import argparse

docker='fanyucai1/virus:latest'
# align reads with bowtie2 and sort bam with samtools
def run(ref_index_dir,outdir,prefix,R1,R2=None):
    outdir=os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)
    ref_index=""
    for i in os.listdir(ref_index_dir):
        if i.endswith(".rev.2.bt2"):
            ref_index = i.split(".rev.2.bt2")[0]
    R1=os.path.abspath(R1)
    cmd=(f'docker run --rm -v {outdir}:/outdir '
         f'-v {R1}:/raw_data/{R1.split("/")[-1]} '
         f'-v {os.path.abspath(ref_index_dir)}:/ref/ ')
    if R2:
        R2=os.path.abspath(R2)
        cmd+=f'-v {R2}:/raw_data/{R2.split("/")[-1]} '
    cmd+=f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bowtie2 --threads 48 -x /ref/{ref_index} '

    if R2:
        cmd+=f'-1 /raw_data/{R1.split("/")[-1]} -2 /raw_data/{R2.split("/")[-1]}|samtools view -bh |samtools sort > /outdir/{prefix}.bam && samtools index /outdir/{prefix}.bam\''
    else:
        cmd= f'-U {R1}|samtools view -bh |samtools sort > /outdir/{prefix}.bam && samtools index /outdir/{prefix}.bam\''
    subprocess.check_call(cmd, shell=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser("Mapping reference.")
    parser.add_argument("-1","--R1",help="R1 fastq",required=True)
    parser.add_argument("-2","--R2",help="R2 fastq",default=None)
    parser.add_argument("-r","--ref",help="directory of bowtie2 reference index",required=True)
    parser.add_argument("-o","--outdir",help="output directory",default=os.getcwd())
    parser.add_argument("-p","--prefix",help="prefix of output")
    args=parser.parse_args()
    run(args.ref,args.outdir,args.prefix,args.R1,args.R2)
