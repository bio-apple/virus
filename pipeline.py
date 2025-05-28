import os
import argparse
import subprocess
import core

docker='virus:latest'
# 获取当前脚本的绝对路径
script_path = os.path.abspath(__file__)

# 获取脚本所在目录
script_dir = os.path.dirname(script_path)

parser=argparse.ArgumentParser("Virus NGS pipeline\nEamile:fanyucai3@gmail.com")
parser.add_argument("-p1","--pe1",help="R1 fastq",required=True, nargs='+')
parser.add_argument("-p2","--pe2",help="R2 fastq",default=None,nargs='+')
parser.add_argument("-p","--prefix",help="prefix of output",required=True, nargs='+')
parser.add_argument("-b","--bed",help="bed file",default=None)
parser.add_argument('-k','--kraken2',help='kraken2 reference index',required=True)
parser.add_argument('-r','--ref',help="directory reference bowtie2 index",required=True)
parser.add_argument('-f','--fna',help="fasta reference",required=True)
parser.add_argument("-h","--host",help="directory host bowtie2 index")
parser.add_argument("-o","--outdir",help="diretory of output",required=True)
args=parser.parse_args()

args.outdir=os.path.abspath(args.outdir)
subprocess.check_call(f'mkdir -p {args.outdir}',shell=True)

cmd=f'docker run --rm -v {script_dir}/modules:/script -v {args.outdir}:/outdir '


for r1,r2,prefix in zip(args.pe1,args.pe2,args.prefix):
    # ------------------------
    # Step 0: fastp qc
    # ------------------------
    r1=os.path.abspath(r1)
    fastp=cmd+f'-v {r1}:/raw_data/{r1.split("/")[-1]} '
    subprocess.check_call(f'mkdir -p {args.outdir}/1.fastp',shell=True)
    if not r2 is None:
        r2 = os.path.abspath(r2)
        fastp += f'-v {r2}:/raw_data/{r2.split("/")[-1]} '
    fastp+=(f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && '
            f'python3 script/fastp.py -R1 /raw_data/{r1.split("/")[-1]} -p {prefix} -o /outdir/1.fastp/ ')
    if not r2 is None:
        fastp+=f'-R2 /raw_data/{r2.split("/")[-1]}'
    fastp+="\'"
    print(fastp)
    subprocess.check_call(cmd,shell=True)

    # ------------------------
    # Step 1: kraken2
    # ------------------------
    args.kraken2=os.path.abspath(args.kraken2)
    kraken2=cmd+(f'-v {os.path.dirname(args.kraken2)}:/ref {docker} sh -c\'export PATH=/opt/conda/envs/kraken2/bin/:$PATH && '
                 f'python3 script/kraken2.py -p1 /outdir/1.fastp/{prefix}.clean_R1.fastq -i /ref/ -o /outdir/ -p {prefix} ')
    if not r2 is None:
        kraken2+=f'-p2 /outdir/1.fastp/{prefix}.clean_R2.fastq '
    kraken2+='\''
    print(kraken2)
    subprocess.check_call(kraken2,shell=True)

    # ------------------------
    # Step 2: bowtie2 host filter
    # ------------------------
    host=os.path.abspath(args.host)
    subprocess.check_call(f'mkdir -p {args.outdir}/2.filter_host', shell=True)
    bowtie2=cmd+(f'-v {host}:/ref/{host.split("/")[-1]} {docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && '
                 f'python3 script/filter_host.py '
                 f'-p1 /outdir/1.fastp/{prefix}.clean_R1.fastq '
                 f'-i /ref/{host.split("/")[-1]} '
                 f'-o /outdir/2.filter_host -p {prefix}')

    if not r2 is None:
        bowtie2+=f'-p2 /outdir/1.fastp/{prefix}.clean_R2.fastq'
    bowtie2+="\'"
    print(bowtie2)
    subprocess.check_call(bowtie2,shell=True)

    # ------------------------
    # Step 3: bowtie2 mapping reference
    # ------------------------
    subprocess.check_call(f'mkdir -p {args.outdir}/3.mapping', shell=True)
    ref=os.path.abspath(args.ref)
    mapping=cmd+(f'-v {ref}:/ref/{ref.split("/")[-1]} {docker} sh -c\'export PATH=/opt/conda/bin/:$PATH && '
                 f'python3 script/fastq2bam.py '
                 f'-1 /outdir/2.filter_host/{prefix}_1.fastq '
                 f'-r /ref/{ref.split("/")[-1]} '
                 f'-o /outdir/3.mapping -p {prefix} ')
    if not r2 is None:
        mapping+=f'-2 /outdir/2.filter_host/{prefix}_2.fastq'
    mapping+="\'"
    print(mapping)
    subprocess.check_call(mapping,shell=True)

    # ------------------------
    # Step 4: trim primer
    # ------------------------
    subprocess.check_call(f'mkdir -p {args.outdir}/4.trim_primer', shell=True)
    if args.bed is not None:
        bed = os.path.abspath(args.bed)
        trim_primer=cmd+(f"{bed}:/ref/{bed.split("/")[-1]} {docker} sh -c\'export PATH=/opt/conda/bin/:$PATH && "
                         f"python3 script/trim_primer -m /outdir/3.mapping/{prefix}.bam "
                         f"-b /ref/{bed.split("/")[-1]} "
                         f"-o /outdir/4.trim_primer/ -p {prefix}\'")
        print(trim_primer)
        subprocess.check_call(trim_primer,shell=True)
    else:
        subprocess.check_call(f'cp {args.outdir}/3.mapping/{prefix}.bam {args.outdir}/4.trim_primer/{prefix}.trimmed.bam',shell=True)

    # ------------------------
    # Step 5: variant calling and consensus sequence
    # ------------------------
    subprocess.check_call(f'mkdir -p {args.outdir}/5.consensus', shell=True)
    args.fna=os.path.abspath(args.fna)
    consensus=cmd+(f"{args.fna}:/ref/{args.fna.split("/")[-1]} {docker} sh -c\'export PATH=/opt/conda/bin/:$PATH && "
                   f"python3 script/consensus.py "
                   f"-b /outdir/4.trim_primer/{prefix}.trimmed.bam "
                   f"-r /ref/{args.fna.split("/")[-1]} "
                   f"-o /outdir/5.consensus -p {prefix}\'")
    print(consensus)
    subprocess.check_call(consensus,shell=True)

    # ------------------------
    # Step 6: plot coverage
    # ------------------------
    subprocess.check_call(f'mkdir -p {args.outdir}/6.coverage',shell=True)
    core.coverage_plot.run(f'{args.outdir}/4.trim_primer/{prefix}.trimmed.bam', f'{args.outdir}/6.coverage', {args.prefix})
    # ------------------------
    # Step 7: run nextclade and pangolin
    # ------------------------
