# Lu J, Rincon N, Wood D E, et al. Metagenome analysis using the Kraken software suite[J]. Nature protocols, 2022, 17(12): 2815-2839.
# Bush S J, Connor T R, Peto T E A, et al. Evaluation of methods for detecting human reads in microbial sequencing datasets[J]. Microbial genomics, 2020, 6(7): e000393.

import os
import subprocess
import argparse

docker="virus:latest"

def run(pe1,outdir,ref,prefix,pe2=None):
    ############prepare input##########################
    pe1=os.path.abspath(pe1)
    outdir=os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)
    ref=os.path.abspath(ref)
    host_index=""
    for i in os.listdir(ref):
        if i.endswith(".rev.2.bt2"):
            host_index=i.split(".rev.2.bt2")[0]
    print(host_index)
    cmd=(f'docker run --rm -v {pe1}:/raw_data/{pe1.split("/")[-1]} '
         f'-v {ref}:/ref '
         f'-v {outdir}:/outdir ')
    if pe2:
        pe2=os.path.abspath(pe2)
        cmd+=f'-v {pe2}:/raw_data/{pe2.split("/")[-1]} '
    cmd+=f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && bowtie2 --very-sensitive-local --no-sq -p 48 -x /ref/{host_index} '

    if pe2:
        pe2 = os.path.abspath(pe2)
        cmd+=(f'-1 /raw_data/{pe1.split("/")[-1]} -2 /raw_data/{pe2.split("/")[-1]} '
              f'--un-conc /outdir/{prefix}_fastq '
              f'-S /outdir/{prefix}.sam > /outdir/{prefix}.bowtie2.log 2>&1\'')
    else:
        cmd+=f'-U /raw_data/{pe1.split("/")[-1]} --un /outdir/{prefix}.unaligned.fastq -S /outdir/{prefix}.sam > /outdir/{prefix}.bowtie2.log 2>&1\''
    print(cmd)
    subprocess.check_call(cmd,shell=True)
    subprocess.check_call(f'rm -rf {outdir}/*.sam',shell=True)
    if pe2 is not None:
        subprocess.check_call(f'mv {outdir}/{prefix}_fastq.1  {outdir}/{prefix}_1.fastq',shell=True)
        subprocess.check_call(f'mv {outdir}/{prefix}_fastq.2  {outdir}/{prefix}_2.fastq',shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Filter host and phix sequence.")
    parser.add_argument("-p1", "--pe1", help="5' reads", required=True)
    parser.add_argument("-p2", "--pe2", help="3' reads", default=None)
    parser.add_argument("-i", "--index", help="directory contains bowtie2 index", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.pe1,args.outdir,args.index,args.prefix,args.pe2)