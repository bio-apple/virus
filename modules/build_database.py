import os,sys,re
import subprocess
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",help="directory of output")
args=parser.parse_args()
args.outdir=os.path.abspath(args.outdir)

docker='fanyucai1/virus:latest'
accession=['NC_004162','NC_045512','NC_001477','NC_001474','NC_001475','NC_002640','NC_063383']
name=['Chikungunya_virus','SARS-CoV-2','Dengue_virus_type_1','Dengue_virus_type_2','Dengue_virus_type_3','Dengue_virus_type_4','Monkeypox_virus']


for a,b in zip(accession,name):
    subprocess.check_call(f'mkdir -p {args.outdir}/{b}',shell=True)
    subprocess.check_call(f'docker run --rm -v {args.outdir}/{b}:/ref/ {docker} sh -c \'export PATH=/opt/conda/bin:$PATH && '
                          f'cd /ref/ && efetch -db nucleotide -id {a} -format fasta >{a}.fasta && '
                          f'bowtie2-build {a}.fasta {a}.fasta && '
                          f'samtools faidx {a}.fasta && '
                          f'bwa index -a bwtsw {a}.fasta\'',shell=True)
