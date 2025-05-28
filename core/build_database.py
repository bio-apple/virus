import os
import subprocess
import argparse

parse=argparse.ArgumentParser()
parse.add_argument("-o","--outdir",help="directory of output",required=True)
args=parse.parse_args()
args.outdir=os.path.abspath(args.outdir)

docker='virus:latest'
accession=['NC_004162','NC_045512','NC_001477','NC_001474','NC_001475','NC_002640','NC_063383']
name=['Chikungunya_virus','SARS-CoV-2','Dengue_virus_type_1','Dengue_virus_type_2','Dengue_virus_type_3','Dengue_virus_type_4','Monkeypox_virus']

# https://genomes.atcc.org/genomes/4f980dee15b2432f
# https://www.ncbi.nlm.nih.gov/nuccore/KX087101
# Zika virus strain ZIKV/Homo sapiens/PRI/PRVABC59/2015, complete genome
name.append("Zika_virus")
accession.append("KX087101")



for a,b in zip(accession,name):
    subprocess.check_call(f'mkdir -p {args.outdir}/{b}',shell=True)
    subprocess.check_call(f'docker run --rm -v {args.outdir}/{b}:/ref/ {docker} sh -c \'export PATH=/opt/conda/envs/kraken2/bin:/opt/conda/bin:$PATH && '
                          f'cd /ref/ && efetch -db nucleotide -id {a} -format fasta >{a}.fasta && '
                          f'bowtie2-build {a}.fasta {a}.fasta && '
                          f'samtools faidx {a}.fasta && '
                          f'bwa index -a bwtsw {a}.fasta\'',shell=True)

