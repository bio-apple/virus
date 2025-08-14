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

parser = argparse.ArgumentParser("Virus NGS pipeline.\nEmail:fanyucai3@gmail.com\n")
# 定义所有共享参数
parser.add_argument("-p1", "--pe1", help="R1 fastq", required=True, nargs='+')
parser.add_argument("-p2", "--pe2", help="R2 fastq", default=None, nargs='+')
parser.add_argument("-p", "--prefix", help="prefix of output", required=True, nargs='+')
parser.add_argument("-o", "--outdir", help="diretory of output", required=True)
parser.add_argument("-c", "--config", help="config file", required=True)
parser.add_argument('-l', '--length', help="read length", type=int, required=True,choices=[50, 75, 100, 150, 200, 250, 300])

# 创建一个参数组，用于组织 ref 和 bowtie2 参数
ref_group = parser.add_argument_group("Reference/Bowtie2 Index Options")
ref_group.add_argument("-e", "--bed", help="bed file", default=None)
ref_group.add_argument("-r", "--ref", help="ref fasta", default=None)
ref_group.add_argument("-b", "--bowtie2", help="directory contains reference bowtie2 index", default=None)
args = parser.parse_args()

# 检查依赖关系
# 依赖关系一：如果 ref 和 bowtie2 只出现一个，则报错
if (args.ref and not args.bowtie2) or (not args.ref and args.bowtie2):
    parser.error("--ref and --bowtie2 must be provided together.")

# 依赖关系二：如果 bed 存在，那么 ref 和 bowtie2 也必须同时存在
if args.bed:
    if not (args.ref and args.bowtie2):
        parser.error("--bed requires both --ref and --bowtie2 to be provided.")
args.outdir=os.path.abspath(args.outdir)
os.makedirs(args.outdir,exist_ok=True)


config = Myconf()
config.read(args.config)

nt_viruses=config.get('database','nt_viruses')
vsp=config.get('database','vsp')
kraken2=config.get('database','kraken2')
host=config.get('database','host')
identify=config.get('parameter','identify')
contig=config.get('parameter','contig_min_length')

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

    read1, read2 = "", ""
    if r2:
        read1 = args.outdir + "/" + "3.filter_host/" + prefix + "_1.fastq"
        read2 = args.outdir + "/" + "3.filter_host/" + prefix + "_2.fastq"
    else:
        read1 = args.outdir + "/" + "3.filter_host/" + prefix + ".unaligned.fastq"
        read2 = None

    if args.ref and args.bowtie2:
        # ------------------------
        # Step 4: mapping
        # ------------------------
        print("#------------------------\n#Step 4: mapping refence\n#------------------------\n")
        core.mapping.run(args.bowtie2,f'{args.outdir}/4.mapping',prefix,read1,read2)

        # ------------------------
        # Step5:trim primer,variant calling,consensus sequence and plot coverage
        # ------------------------
        print("#------------------------\n#Step5:trim primer,variant calling,consensus sequence and plot coverage\n#------------------------\n")
        if args.bed:
            # trim primer
            core.trim_primer.run(args.bed, f'{args.outdir}/4.mapping/{prefix}.bam', f'{args.outdir}/5.consensus/ref/',
                                 prefix)
            # consensus
            core.consensus.run(f'{args.outdir}/5.consensus/ref/{prefix}.soft.clipped.sort.bam', f'{args.outdir}/5.consensus/ref/', prefix,args.ref)
        else:
            # consensus
            core.consensus.run(f'{args.outdir}/4.mapping/{prefix}.bam',f'{args.outdir}/5.consensus/ref/', prefix,args.ref)
    else:
        # ------------------------
        # Step 4: denovo genome assembly(megahit and metaspades) and remove redundancy (cd-hit-est)
        # ------------------------
        print("#------------------------\n#Step 4: denovo genome assembly(megahit and metaspades) and remove redundancy (cd-hit-est)\n#------------------------\n")
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [
                executor.submit(core.megahit.run, read1, prefix, args.outdir + "/4.assembly/", read2,contig),
                executor.submit(core.metaspades.run, read1, prefix, args.outdir + "/4.assembly/", read2)
            ]
            for future in as_completed(futures):
                print(future.result())
        subprocess.check_call(f'cd {args.outdir}/4.assembly/ && cat spades_{prefix}/scaffolds_{contig}bp.fasta megahit_{prefix}/{prefix}.contigs.fa >{prefix}.contigs.fa',shell=True)
        core.cd_hit_est.run(f'{args.outdir}/4.assembly/{prefix}.contigs.fa',identify,prefix+".non-redundant",f'{args.outdir}/4.assembly/')
        #build bowtie2 index
        index = (f'docker run --rm -v {args.outdir}/4.assembly/:/raw_data/ {docker} sh -c '
                 f'\'export PATH=/opt/conda/bin/:$PATH && '
                 f'bowtie2-build /raw_data/{prefix}.non-redundant.fna /raw_data/{prefix}.non-redundant.fna\' ')

        print(index)
        subprocess.check_call(index, shell=True)

        # ------------------------
        # Step 5: blast NCBI Database: nt virus and parse blast result
        # ------------------------
        print("#------------------------\n#Step 5: blast NCBI Database: nt virus and vsp\n#------------------------\n")
        core.blast.run(f'{args.outdir}/4.assembly/{prefix}.non-redundant.fna',nt_viruses,f"{args.outdir}/5.blast/",prefix+".nt_viruses",10)
        core.blast.run(f'{args.outdir}/4.assembly/{prefix}.non-redundant.fna', vsp, f"{args.outdir}/5.blast/", prefix+".vsp",10)
        num=core.parse_blast.run(f"{args.outdir}/5.blast/{prefix}.vsp.blast_all.txt",f"{args.outdir}/5.blast/{prefix}.nt_viruses.blast_all.txt",nt_viruses,f"{args.outdir}/5.blast/")
        # ------------------------
        # step 6:mapping reference
        # ------------------------
        print("#------------------------\n#Step 6:mapping reference\n#------------------------\n")
        if num!=0:
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
        core.consensus.run(f'{args.outdir}/6.mapping/denovo/{prefix}.bam', f'{args.outdir}/7.consensus/denovo', prefix)
        if num!=0:
            core.consensus.run(f'{args.outdir}/7.consensus/ref/{prefix}.soft.clipped.sort.bam',f'{args.outdir}/7.consensus/ref/', prefix)