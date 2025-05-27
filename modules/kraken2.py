#Lu J, Rincon N, Wood D E, et al. Metagenome analysis using the Kraken software suite[J]. Nature protocols, 2022, 17(12): 2815-2839.
#Shen Z, Robert L, Stolpman M, et al. A genome catalog of the early-life human skin microbiome[J]. Genome Biology, 2023, 24(1): 252.
import sys
import os
import subprocess
import argparse
import gzip

def run(pe1,ref,prefix,outdir,pe2=None):
    pe1 = os.path.abspath(pe1)
    outdir = os.path.abspath(outdir)
    read_length=get_average_read_length(pe1)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cmd = f"kraken2 --confidence 0.6 --db {ref} --threads 24 --output {outdir}/{prefix}.txt --minimum-base-quality 20 --report {outdir}/{prefix}.report.txt "
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

def get_average_read_length(fastq_file, max_reads=10):
    # 判断是否是压缩文件
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    total_len = 0
    read_count = 0

    with open_func(fastq_file, "rt") as f:
        while read_count < max_reads:
            header = f.readline()
            if not header:
                break  # 文件结束
            seq = f.readline().strip()     # 序列
            f.readline()                   # 跳过 '+'
            f.readline()                   # 跳过质量值

            total_len += len(seq)
            read_count += 1

    if read_count == 0:
        return 0
    read_length=total_len / read_count
    if  read_length in(301,151,101,251,51,201):
        read_length=read_length-1
    return int(read_length)  # 返回平均长度


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Classified out option on the Kraken database,")
    parser.add_argument("-p1", "--pe1", help="5' reads", required=True)
    parser.add_argument("-p2", "--pe2", help="3' reads",default=None)
    parser.add_argument("-i", "--index", help="directory contains kraken2 index", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.pe1,args.index,args.prefix,args.outdir,args.pe2)
