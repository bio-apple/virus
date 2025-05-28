#   Lu J, Rincon N, Wood D E, et al. Metagenome analysis using the Kraken software suite[J]. Nature protocols, 2022, 17(12): 2815-2839.
#   Shen Z, Robert L, Stolpman M, et al. A genome catalog of the early-life human skin microbiome[J]. Genome Biology, 2023, 24(1): 252.
#   Liu Y, Ghaffari M H, Ma T, et al. Impact of database choice and confidence score on the performance of taxonomic classification using Kraken2[J]. aBIOTECH, 2024: 1-11.

import os,gzip
import subprocess
import argparse

docker="virus:latest"
def run(pe1,index,prefix,outdir,pe2=None):
    pe1 = os.path.abspath(pe1)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    cmd=f"docker run --rm -v {pe1}:/raw_data/{os.path.basename(pe1)} -v {os.path.abspath(index)}:/ref -v {outdir}:/outdir "

    if pe2:
        pe2=os.path.abspath(pe2)
        cmd+=f"-v {pe2}:/raw_data/{os.path.basename(pe2)} "

    read_length=get_average_read_length(pe1)
    cmd+=f"{docker} sh -c \'export PATH=/opt/conda/envs/kraken2/bin/:$PATH && kraken2 --confidence 0.4 --db /ref --threads 24 --output /outdir/{prefix}.txt --minimum-base-quality 20 --report /outdir/{prefix}.report.txt "
    if pe2:
        cmd+=f"{'--paired /raw_data/' + os.path.basename(pe1) + ' /raw_data/' + os.path.basename(pe2) if pe2 else '/raw_data/' + os.path.basename(pe1)} && "

    # Run Bracken for Abundance Estimation of Microbiome Samples
    cmd += f"bracken -d /ref/ -i /outdir/{prefix}.report.txt -r {read_length} -o /outdir/{prefix}.bracken -w /outdir/{prefix}.breport -t 10 && "

    # Generate Krona Plots
    cmd += f"kreport2krona.py -r /outdir/{prefix}.breport -o /outdir/{prefix}.krona.txt --no-intermediate-ranks && ktImportText /outdir/{prefix}.krona.txt -o /outdir/{prefix}.krona.html\'"
    print(cmd)
    subprocess.check_call(cmd,shell=True)

def get_average_read_length(fastq_file, max_reads=10):
    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, "rt") as f:
        lengths = [len(f.readline().strip()) for _ in range(max_reads)
                   if f.readline()]  # 跳 header，读取序列后跳过+和质量
        [f.readline() for _ in range(2 * len(lengths))]

    if not lengths:
        return 0

    avg = sum(lengths) / len(lengths)
    return int(avg - 1) if int(avg) in {301, 151, 101, 251, 51, 201} else int(avg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Classified out option on the Kraken database,")
    parser.add_argument("-p1", "--pe1", help="5' reads", required=True)
    parser.add_argument("-p2", "--pe2", help="3' reads",default=None)
    parser.add_argument("-i", "--index", help="directory contains kraken2 index", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.pe1,args.index,args.prefix,args.outdir,args.pe2)
