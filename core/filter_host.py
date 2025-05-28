# Lu J, Rincon N, Wood D E, et al. Metagenome analysis using the Kraken software suite[J]. Nature protocols, 2022, 17(12): 2815-2839.
# Bush S J, Connor T R, Peto T E A, et al. Evaluation of methods for detecting human reads in microbial sequencing datasets[J]. Microbial genomics, 2020, 6(7): e000393.

import os
import subprocess
import argparse

docker="virus:latest"
parser=argparse.ArgumentParser("Filter human host and phix sequence.")
parser.add_argument("-p1","--pe1",help="5' reads",required=True)
parser.add_argument("-p2","--pe2",help="3' reads",default=None)
parser.add_argument("-i","--index",help="directory contains bowtie2 index",required=True)
parser.add_argument("-o","--outdir",help="output directory",required=True)
parser.add_argument("-p","--prefix",help="prefix of output",required=True)
args=parser.parse_args()

############prepare input##########################
args.pe1=os.path.abspath(args.pe1)
in_dir=os.path.dirname(args.pe1)
a=args.pe1.split("/")[-1]
args.outdir=os.path.abspath(args.outdir)
if not os.path.exists(args.outdir):
    subprocess.check_call('mkdir -p %s'%(args.outdir),shell=True)
args.index=os.path.abspath(args.index)
ref=os.path.abspath(args.index)
host_index=""
for i in os.listdir(args.index):
    if i.endswith(".rev.2.bt2"):
        host_index=i.split(".rev.2.bt2")[0]
print(host_index)
cmd=f"docker run -v {in_dir}:/raw_data/ -v {ref}:/ref -v {args.outdir}:/outdir {docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && "
if args.pe2 is not None:
    args.pe2 = os.path.abspath(args.pe2)
    if in_dir!=os.path.dirname(args.pe2):
        print("read1 and reads2 must be in the same directory.")
        exit()
    b=args.pe2.split("/")[-1]
    cmd+=(f"bowtie2 --very-sensitive-local --no-sq -p 48 -x /ref/{host_index} "
         f"-1 /raw_data/{a} -2 /raw_data/{b} "
         f"--un-conc /outdir/{args.prefix}_fastq "
         f"-S /outdir/{args.prefix}.sam > /outdir/{args.prefix}.bowtie2.log 2>&1\'")
else:
    cmd+=f"bowtie2 --very-sensitive-local -p 48 -U /raw_data/{a} -x /ref/{host_index} --un /outdir/{args.prefix}.unaligned.fastq -S /outdir/{args.prefix}.sam > /outdir/{args.prefix}.bowtie2.log 2>&1\'"
print(cmd)
subprocess.check_call(cmd,shell=True)
subprocess.check_call(f'rm -rf {args.outdir}/*.sam',shell=True)
if args.pe2 is not None:
    subprocess.check_call(f'mv {args.outdir}/{args.prefix}_fastq.1  {args.outdir}/{args.prefix}_1.fastq',shell=True)
    subprocess.check_call(f'mv {args.outdir}/{args.prefix}_fastq.2  {args.outdir}/{args.prefix}_2.fastq',shell=True)