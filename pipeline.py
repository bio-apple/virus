import os
import argparse
import subprocess
import configparser
import core
from concurrent.futures import ThreadPoolExecutor, as_completed

class Myconf(configparser.ConfigParser):
    def __init__(self, defaults=None):
        configparser.ConfigParser.__init__(self, defaults=defaults)
    def optionxform(self, optionstr):
        return optionstr

docker='virus:latest'
# 获取当前脚本的绝对路径
script_path = os.path.abspath(__file__)

# 获取脚本所在目录
script_dir = os.path.dirname(script_path)

#get vsp species list
accession=os.path.join(script_dir,'/plugin/accession.list')
ssname=os.path.join(script_dir,'/plugin/sscinames.txt')

parser=argparse.ArgumentParser("Virus NGS pipeline.\nEmail:fanyucai3@gmail.com\n")
parser.add_argument("-p1","--pe1",help="R1 fastq",required=True, nargs='+')
parser.add_argument("-p2","--pe2",help="R2 fastq",default=None,nargs='+')
parser.add_argument("-p","--prefix",help="prefix of output",required=True, nargs='+')
parser.add_argument("-o","--outdir",help="diretory of output",required=True)
parser.add_argument("-c","--config",help="config file",required=True)
parser.add_argument('-l','--length',help="read length",type=int,required=True,choices=[50,75,100,150,200,250,300])
args=parser.parse_args()

args.outdir=os.path.abspath(args.outdir)
os.makedirs(args.outdir,exist_ok=True)


config = Myconf()
config.read(args.config)

blastdb=config.get('database','nt_viruses')
kraken2=config.get('database','kraken2')
host=config.get('database','host')
identify=config.get('parameter','identify')
contig=config.get('parameter','contig_min_length')
bed=config.get('bed_file','bed')

for r1,r2,prefix in zip(args.pe1,args.pe2,args.prefix):
    # ------------------------
    # Step 1: fastp qc
    # ------------------------
    print("#------------------------\n#Step1:fastp qc\n#------------------------\n")
    core.fastp.run(r1,args.outdir+"/1.fastp",prefix,r2)

    # ------------------------
    # Step 2: kraken2
    # ------------------------
    print("#------------------------\n#Step 2: kraken2\n#------------------------\n")
    core.kraken2.run(r1,kraken2,prefix,args.outdir+"/2.kraken2",args.length,r2)

    # ------------------------
    # Step 3: bowtie2 host filter
    # ------------------------
    print("#------------------------\n#Step 3: bowtie2 host filter\n#------------------------\n")
    core.filter_host.run(r1,args.outdir+"/3.filter_host",host,prefix,r2)

    # ------------------------
    # Step 4: denovo genome assembly(megahit and metaspades) and remove redundancy (cd-hit-est)
    # ------------------------
    print("#------------------------\n#Step 4: denovo genome assembly(megahit and metaspades) and remove redundancy (cd-hit-est)\n#------------------------\n")
    read1,read2="",""
    if r2:
        read1=args.outdir+"/"+"3.filter_host/"+prefix+"_1.fastq"
        read2=args.outdir+"/"+"3.filter_host/"+prefix+"_2.fastq"
    else:
        read1 = args.outdir + "/" + "3.filter_host/" + prefix +".unaligned.fastq"
        read2=None
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(core.megahit.run, read1, prefix, args.outdir + "/4.assembly/", read2,contig),
            executor.submit(core.metaspades.run, read1, prefix, args.outdir + "/4.assembly/", read2)
        ]
        for future in as_completed(futures):
            print(future.result())
    subprocess.check_call(f'cd {args.outdir}/4.assembly/ && cat spades_{prefix}/scaffolds_{contig}bp.fasta megahit_{prefix}/{prefix}.contigs.fa >{prefix}.contigs.fa',shell=True)
    core.cd_hit_est.run(f'{args.outdir}/4.assembly/{prefix}.contigs.fa',identify,prefix+".non-redundant",f'{args.outdir}/4.assembly/')

    # ------------------------
    # Step 5: blast NCBI Database: nt virus and parse blast result and find corresponding species in VSPv2
    # ------------------------
    print("#------------------------\n#Step 5: blast NCBI Database: nt virus\n#------------------------\n")
    core.blast.run(f'{args.outdir}/4.assembly/{prefix}.non-redundant.fna',blastdb,f"{args.outdir}/5.blast/",prefix,10)
    index=(f'docker run --rm -v {args.outdir}/4.assembly/:/raw_data/ {docker} sh -c '
           f'\'export PATH=/opt/conda/bin/:$PATH && '
           f'bowtie2-build /raw_data/{prefix}.non-redundant.fna /raw_data/{prefix}.non-redundant.fna\' ')

    print(index)
    subprocess.check_call(index,shell=True)
    core.blast2VSP.run(f"{args.outdir}/5.blast/{prefix}.blast_all.txt",accession,ssname,blastdb,f"{args.outdir}/5.blast/")
    chr = []
    infile = open(f"{args.outdir}/5.blast/{prefix}.blast_all.txt", "r")
    for line in infile:
        line = line.strip()
        if not line.startswith("#"):
            array = line.split("\t")
            if not array[0] in chr:
                chr.append(array[0])
    print(chr)

    # ------------------------
    # step 6:mapping reference
    # ------------------------
    print("#------------------------\n#Step 6:mapping reference\n#------------------------\n")
    if os.path.exists(f"{args.outdir}/5.blast/ref.fasta"):
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [
                executor.submit(core.mapping.run, f'{args.outdir}/5.blast/',f'{args.outdir}/6.mapping/ref',prefix,r1, r2),
                executor.submit(core.mapping.run, f'{args.outdir}/4.assembly/',f'{args.outdir}/6.mapping/denovo',prefix,r1,r2)
            ]
            for future in as_completed(futures):
                print(future.result())
    else:
        core.mapping.run(f'{args.outdir}/4.assembly/',f'{args.outdir}/6.mapping/denovo',prefix,r1,r2)

    # ------------------------
    # step7:trim primer,variant calling,consensus sequence and plot coverage
    # ------------------------
    print("#------------------------\n#Step7:trim primer,variant calling,consensus sequence and plot coverage\n#------------------------\n")
    core.consensus.run(f'{args.outdir}/6.mapping/denovo/{prefix}.bam', f'{args.outdir}/7.consensus/denovo', prefix,None, " ".join(chr))
    if os.path.exists(f"{args.outdir}/5.blast/ref.fasta"):
        if bed:
            core.trim_primer.run(bed,f'{args.outdir}/6.mapping/ref/{prefix}.bam', f'{args.outdir}/7.consensus/ref/',prefix)
            core.consensus.run(f'{args.outdir}/7.consensus/ref/{prefix}.soft.clipped.sort.bam', f'{args.outdir}/7.consensus/ref/', prefix,args.ref)
        else:
            core.consensus.run(f'{args.outdir}/6.mapping/ref/{prefix}.bam', f'{args.outdir}/7.consensus/ref/', prefix, f"{args.outdir}/5.blast/ref.fasta")
