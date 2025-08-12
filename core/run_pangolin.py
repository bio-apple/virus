import os
import subprocess
import argparse

docker="virus:latest"

def run(ref,outdir,prefix):
    ref=os.path.abspath(ref)
    outdir=os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)
    cmd = (f"docker run "
           f"-v {ref}/:/ref/{ref.split("/")[-1]} "
           f"-v {outdir}:/outdir/ {docker} "
           f"sh -c \'export PATH=/opt/conda/envs/pangolin/bin/:$PATH && "
           f"pangolin /ref/{ref.split("/")[-1]} --threads 8 --outfile /outdir/{prefix}.pangolin.results.csv\'")

    print(cmd)
    subprocess.run(cmd, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Run pangolin analysis.")
    parser.add_argument("-r","--ref", help="Reference fasta file",required=True)
    parser.add_argument("-o","--outdir", help="Output directory",default=os.getcwd())
    parser.add_argument("-p","--prefix", help="prefix of output files",required=True)
    args = parser.parse_args()
    run(args.ref,args.outdir,args.prefix)