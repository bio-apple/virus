#Lu J, Rincon N, Wood D E, et al. Metagenome analysis using the Kraken software suite[J]. Nature protocols, 2022, 17(12): 2815-2839.
#Shen Z, Robert L, Stolpman M, et al. A genome catalog of the early-life human skin microbiome[J]. Genome Biology, 2023, 24(1): 252.
import sys
import os
import subprocess
import argparse

def run(pe1,ref,prefix,outdir,pe2=None,read_length=150):
    pe1 = os.path.abspath(pe1)
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cmd = f"kraken2 --confidence 0.4 --db {ref} --threads 24 --output {outdir}/{prefix}.txt --minimum-base-quality 20 --report {outdir}/{prefix}.report.txt "
    if pe2 is not None:
        pe2=os.path.abspath(pe2)
        cmd+=f"--paired {pe1} {pe2}"
    else:
        cmd+=f"{pe1}"
    #Run Bracken for Abundance Estimation of Microbiome Samples
    cmd+=f" && bracken -d {ref} -i {outdir}/{prefix}.report.txt -r {read_length} -o {outdir}/{prefix}.bracken -w {outdir}/{prefix}.breport -t 10"
    #Generate Krona Plots
    cmd+=(f" && kreport2krona.py -r {outdir}/{prefix}.breport -o {outdir}/{prefix}.krona.txt --no-intermediate-ranks && "
          f"ktImportText {outdir}/{prefix}.krona.txt -o {outdir}/{prefix}.krona.html")
    print(cmd)
    subprocess.check_call(cmd,shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Classified out option on the Kraken database,")
    parser.add_argument("-p1", "--pe1", help="5' reads", required=True)
    parser.add_argument("-r","--read_length", help="Read length", required=True,default=150)
    parser.add_argument("-p2", "--pe2", help="3' reads",default=None)
    parser.add_argument("-i", "--index", help="directory contains kraken2 index", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.pe1,args.index,args.prefix,args.outdir,args.pe2,args.read_length)
