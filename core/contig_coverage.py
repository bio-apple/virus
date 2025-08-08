import os,sys,re
import argparse
import subprocess

docker='virus:latest'
def run(ref,pe1,outdir,prefix,pe2=None,cov=5):
    ref=os.path.abspath(ref)
    pe1=os.path.abspath(pe1)
    pe2=os.path.abspath(pe2)
    outdir=os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)
    if pe2 is not None:
        cmd=(f'docker run '
             f'-v {pe1}:/raw_data/{pe1.split("/")[-1]} '
             f'-v {pe2}:/raw_data/{pe2.split("/")[-1]} '
             f'-v {outdir}:/outdir/ '
             f'-v {ref}:/ref/{ref.split("/")[-1]} {docker} '
             f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
             f'bbwrap.sh ref=/ref/{ref.split("/")[-1]} in1=/raw_data/{pe1.split("/")[-1]} in2=/raw_data/{pe2.split("/")[-1]} out=/outdir/{prefix}.sam.gz && '
             f'pileup.sh in=/outdir/{prefix}.sam.gz out=/outdir/{prefix}.cov.txt\'')
    else:
        cmd = (f'docker run '
               f'-v {pe1}:/raw_data/{pe1.split("/")[-1]} '
               f'-v {outdir}:/outdir/ '
               f'-v {ref}:/ref/{ref.split("/")[-1]} '
               f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
               f'bbwrap.sh ref=/ref/{ref.split("/")[-1]} in=/raw_data/{pe1.split("/")[-1]} out=/outdir/{prefix}.sam.gz && '
               f'pileup.sh in=/outdir/{prefix}.sam.gz out=/outdir/{prefix}.cov.txt\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Calculate coverage of contigs')
    parser.add_argument("-r",'--ref',help='Reference fasta file',required=True)
    parser.add_argument("-p1","--pe1",help='PE1 fastq file',required=True)
    parser.add_argument("-p2","--pe2",help='PE2 fastq file',default=None)
    parser.add_argument("-o","--outdir",help='Output directory',default=os.getcwd())
    parser.add_argument("-p","--prefix",help='Output prefix',required=True)
    parser.add_argument("-c","--coverage",help='min coverage threshold,default=5',default=5,type=int,required=True)
    args = parser.parse_args()
    run(args.ref,args.pe1,args.outdir,args.prefix,args.pe2,args.coverage)