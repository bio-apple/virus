import os,sys,re
import argparse
import subprocess

docker='virus:latest'

# 获取当前脚本的绝对路径
script_path = os.path.abspath(__file__)

# 获取脚本所在目录
script_dir = os.path.dirname(script_path)

parser=argparse.ArgumentParser("Virus NGS pipeline\nEamile:fanyucai3@gmail.com")
parser.add_argument("-p1","--pe1",help="R1 fastq",required=True, nargs='+')
parser.add_argument("-p2","--pe2",help="R2 fastq",default=None,nargs='+')
parser.add_argument("-p","--prefix",help="prefix of output",required=True, nargs='+')
parser.add_argument("-b","--bed",help="bed file")
parser.add_argument("-o","--outdir",help="diretory of output",required=True)
args=parser.parse_args()

args.outdir=os.path.abspath(args.outdir)
subprocess.check_call(f'mkdir -p {args.outdir}',shell=True)

cmd=f'docker run --rm -v {script_dir}/modules:/script -v {args.outdir}:/outdir '



for a,b,c in zip(args.pe1,args.pe2,args.prefix):
    a=os.path.abspath(a)
    if not b is None:
        b=os.path.abspath(a)
    #step1:fastp
    step1=cmd+(f'-v {a}:/raw_data/{a.split("/")[-1]} -v {b}:/raw_data/{b.split("/")[-1]} {docker} sh -c \''
               f'export PATH=/opt/conda/bin/:$PATH && '
               f'python3 script/fastp.py -R1 /raw_data/{a.split("/")[-1]} '
               f'-R2 /raw_data/{b.split("/")[-1]} '
               f'-p {c} -o /outdir/1.fastp/ \'')
    subprocess.check_call(cmd,shell=True)

    #step2:filter host
    step2=cmd+(f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && '
               f'python3 script/filter_host.py "\'')


