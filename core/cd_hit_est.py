# Eamil:yucai.fan@illumina.com
import os
import subprocess
import argparse
import time

docker="virus:latest"
def run(fna,identify,prefix,outdir):
    start=time.time()
    outdir=os.path.abspath(outdir)
    fna=os.path.abspath(fna)
    os.makedirs(outdir,exist_ok=True)
    # Lythgoe K A, Hall M, Ferretti L, et al. SARS-CoV-2 within-host diversity and transmission[J]. Science, 2021, 372(6539): eabg0821.
    # 0.995 sequence identity
    # Dezordi F Z, Resende P C, Naveca F G, et al. Unusual SARS-CoV-2 intrahost diversity reveals lineage superinfection[J]. Microbial Genomics, 2022, 8(3).
    # 99.8% sequence identity
    cmd=(f'docker run --rm -v {fna}:/raw_data/{fna.split("/")[-1]} -v {outdir}:/outdir {docker} sh -c \''
         f'export PATH=/opt/conda/bin:$PATH && '
         f'cd-hit-est -i /raw_data/{fna.split("/")[-1]} -o /outdir/{prefix}.fna -c {identify}\'')
    subprocess.check_call(cmd,shell=True)

    end=time.time()
    print("Elapse time is %g seconds" % (end - start))

if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("-f", "--fna", help="fasta sequence", required=True)
    parser.add_argument("-c", "--identify", help="sequence identity threshold, default: 0.998", default=0.998,
                        type=float, required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", default=time.strftime("%Y-%m-%d"))
    parser.add_argument("-o", "--outdir", help="output directory", default=os.getcwd(), required=True)
    args = parser.parse_args()
    run(args.fna,args.identify,args.prefix,args.outdir)