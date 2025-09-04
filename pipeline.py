import os
import argparse
import subprocess
import configparser
import core
from concurrent.futures import ThreadPoolExecutor, as_completed
import time


class Myconf(configparser.ConfigParser):
    def __init__(self, defaults=None):
        configparser.ConfigParser.__init__(self, defaults=defaults)
    def optionxform(self, optionstr):
        return optionstr

docker='virus:latest'
# Get the absolute path of the current script
script_path = os.path.abspath(__file__)

# Get the script's directory
script_dir = os.path.dirname(script_path)

parser = argparse.ArgumentParser("Virus NGS pipeline.\nEmail:fanyucai3@gmail.com\n")
# Define all shared parameters
parser.add_argument("-p1", "--pe1", help="R1 fastq", required=True, nargs='+')
parser.add_argument("-p2", "--pe2", help="R2 fastq", default=None, nargs='+')
parser.add_argument("-p", "--prefix", help="prefix of output", required=True, nargs='+')
parser.add_argument("-o", "--outdir", help="diretory of output", required=True)
parser.add_argument("-c", "--config", help="config file", required=True)
parser.add_argument('-l', '--length', help="read length", type=int, required=True,choices=[50, 75, 100, 150, 200, 250, 300])

# Create a parameter group to organize the ref and bowtie2 parameters.
ref_group = parser.add_argument_group("Reference/Bowtie2 Index Options")
ref_group.add_argument("-e", "--bed", help="bed file(Optional)", default=None)
ref_group.add_argument("-r", "--ref", help="ref fasta(Optional)", default=None)
ref_group.add_argument("-b", "--bowtie2", help="directory contains reference bowtie2 index(Optional)", default=None)
args = parser.parse_args()

##########################################################################Dependency Check
#Dependency 1: If only one of 'ref' and 'bowtie2' exists, report an error.
if (args.ref and not args.bowtie2) or (not args.ref and args.bowtie2):
    parser.error("--ref and --bowtie2 must be provided together.")

#Dependency 2: If bed exists, then ref and bowtie2 must also exist simultaneously.
if args.bed:
    if not (args.ref and args.bowtie2):
        parser.error("--bed requires both --ref and --bowtie2 to be provided.")
args.outdir=os.path.abspath(args.outdir)
os.makedirs(args.outdir,exist_ok=True)


config = Myconf()
config.read(args.config)

virus=config.get('database','virus')
nt_viruses=config.get('blast_db','nt_viruses')
vsp=config.get('blast_db','vsp')
kraken2=config.get('database','kraken2')
host=config.get('database','host')
vsp_fa=config.get('fasta','vsp')
identify=config.get('parameter','identify')
contig=config.get('parameter','contig_min_length')

for r1,r2,prefix in zip(args.pe1,args.pe2,args.prefix):
    start = time.time()

    # ------------------------
    # Step 1: fastp qc
    # ------------------------
    print("\n#------------------------\n#Step1:fastp qc\n#------------------------\n")
    core.fastp.run(r1,args.outdir+"/1.fastp",prefix,r2)

    # ------------------------
    # Step 2: kraken2
    # ------------------------
    print("\n#------------------------\n#Step 2: kraken2\n#------------------------\n")
    core.kraken2.run(r1,kraken2,prefix,args.outdir+"/2.kraken2",args.length,r2)

    # ------------------------
    # Step 3: bowtie2 host filter
    # ------------------------
    print("\n#------------------------\n#Step 3: bowtie2 host filter\n#------------------------\n")
    core.filter_host.run(r1,args.outdir+"/3.filter_host",host,prefix,r2)

    read1, read2 ,accession1,count= "", None,[],0
    if r2 is None:
        count=core.subsample.run(args.outdir + "/" + "3.filter_host/" + prefix + ".unaligned.fastq.gz",args.outdir+"/3.filter_host",prefix)
        if count > 4000000:
            read1 = args.outdir + "/" + "3.filter_host/" + prefix + ".sub.R1.fastq.gz"
        else:
            read1 = args.outdir + "/" + "3.filter_host/" + prefix + ".unaligned.fastq.gz"
    else:
        read1 = args.outdir + "/" + "3.filter_host/" + prefix + "_1.fastq.gz"
        read2 = args.outdir + "/" + "3.filter_host/" + prefix + "_2.fastq.gz"
        count = core.subsample.run(read1,args.outdir + "/3.filter_host", prefix,read2)
        if count>4000000:
            read1 = args.outdir + "/" + "3.filter_host/" + prefix + ".sub.R1.fastq.gz"
            read2 = args.outdir + "/" + "3.filter_host/" + prefix + ".sub.R2.fastq.gz"
        if args.length >100 and not args.ref:
            count=core.merge_fastq.run(read1, read2, prefix, f'{args.outdir}/3.filter_host/', args.length)
            accession1 = core.blast.run(f'{args.outdir}/3.filter_host/{prefix}.non-redundant.fna', virus,f'{args.outdir}/3.filter_host/', prefix, 10, 98, 95, 1e-10, 1)
            print(accession1)
    if args.ref and args.bowtie2:
        os.makedirs(f"{args.outdir}/4.vsp",exist_ok=True)
        os.makedirs(f"{args.outdir}/5.assembly/", exist_ok=True)
        os.makedirs(f"{args.outdir}/6.blast/", exist_ok=True)

        # ------------------------
        # Step 7: mapping
        # ------------------------
        print("\n#------------------------\n#Step 7: mapping refence\n#------------------------\n")
        core.mapping.run(args.bowtie2,f'{args.outdir}/7.mapping',prefix,read1,read2)

        # ------------------------
        # Step8:trim primer,variant calling,consensus sequence and plot coverage
        # ------------------------
        print("#------------------------\n#Step8:trim primer,variant calling,consensus sequence and plot coverage\n#------------------------\n")
        if args.bed:
            # trim primer
            core.trim_primer.run(args.bed, f'{args.outdir}/7.mapping/{prefix}.bam', f'{args.outdir}/8.consensus/',prefix)
            # consensus
            core.ref_consensus.run(f'{args.outdir}/8.consensus/{prefix}.soft.clipped.sort.bam', f'{args.outdir}/8.consensus/', prefix,args.ref)
        else:
            # consensus
            core.ref_consensus.run(f'{args.outdir}/7.mapping/{prefix}.bam',f'{args.outdir}/8.consensus/', prefix,args.ref)
    else:

        # ------------------------
        # Step 4: non-host reads mapping VSP database
        # ------------------------
        print("\n#------------------------\n#Step 4: non-host reads mapping VSP database\n#------------------------\n")
        accession2=core.contig_cov.run(vsp_fa,read1,f'{args.outdir}/4.vsp/',prefix,read2,args.length)
        accession = list(set(accession1 + accession2))
        print(accession)

        # ------------------------
        # Step 5: denovo genome assembly(megahit and metaspades) and remove redundancy (cd-hit-est)
        # ------------------------
        print("\n#------------------------\n#Step 5: denovo genome assembly(megahit and metaspades) and remove redundancy (cd-hit-est)\n#------------------------\n")
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [
                executor.submit(core.megahit.run, read1, prefix, args.outdir + "/5.assembly/", read2,contig),
                executor.submit(core.metaspades.run, read1, prefix, args.outdir + "/5.assembly/", read2)
            ]
            for future in as_completed(futures):
                print(future.result())
        subprocess.check_call(f'cd {args.outdir}/5.assembly/ && cat spades_{prefix}/scaffolds_{contig}bp.fasta megahit_{prefix}/{prefix}.contigs.fa >{prefix}.contigs.fa',shell=True)
        core.cd_hit_est.run(f'{args.outdir}/5.assembly/{prefix}.contigs.fa',identify,prefix+".non-redundant",f'{args.outdir}/5.assembly/')

        # ------------------------
        # Step 6: blast nt
        # ------------------------
        print("\n#------------------------\n#Step 6: blast NCBI Database: nt virus\n#------------------------\n")
        core.blast.run(f'{args.outdir}/5.assembly/{prefix}.non-redundant.fna',virus,f"{args.outdir}/6.blast/",prefix+".nt_viruses",10,98,70,1e-10,1)
        new_accession1 = core.parse_blast.run(f"{args.outdir}/6.blast/{prefix}.nt_viruses.blast_all.txt",nt_viruses, f"{args.outdir}/6.blast/", accession,0.95)
        print(new_accession1)

        # ------------------------
        # step 7:mapping reference
        # ------------------------
        print("\n#------------------------\n#Step 7:mapping reference\n#------------------------\n")
        final_accession=core.mapping.run(f'{args.outdir}/6.blast/',f'{args.outdir}/7.mapping/',prefix,read1, read2,args.length)
        print(final_accession)

        # ------------------------
        # step8:variant calling,consensus sequence and plot coverage
        # ------------------------
        print("#------------------------\n#Step7:variant calling,consensus sequence and plot coverage\n#------------------------\n")
        core.consensus.run(final_accession,read1,f'{args.outdir}/8.consensus/',prefix,nt_viruses,args.length,read2)

    end=time.time()
    print(f"\nSampleID {prefix}:Elapse time is {(end-start)} seconds\n")