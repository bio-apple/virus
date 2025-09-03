import os,sys
import subprocess,argparse

docker="virus:latest"

def run(pe1,pe2,prefix,outdir,length=150):
    pe1=os.path.abspath(pe1)
    pe2=os.path.abspath(pe2)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)
    #assmebly
    cmd=(f'docker run --rm -v {outdir}:/outdir '
         f'-v {pe1}:/raw_data/{pe1.split("/")[-1]} '
         f'-v {pe2}:/raw_data/{pe2.split("/")[-1]} '
         f'{docker} sh -c \'export PATH=/opt/conda/bin/:$PATH && '
         f'bbmerge.sh in1=/raw_data/{pe1.split("/")[-1]} in2=/raw_data/{pe2.split("/")[-1]} out=/outdir/{prefix}.tmp.fasta && '
         f'reformat.sh in=/outdir/{prefix}.tmp.fasta out=/outdir/{prefix}.merge.fasta minlength={length}\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    #cd-hit-est
    cmd = (f'docker run --rm -v {outdir}:/outdir {docker} sh -c \''
           f'export PATH=/opt/conda/bin:$PATH && '
           f'cd-hit-est -i /outdir/{prefix}.merge.fasta -o /outdir/{prefix}.non-redundant.fna -M 0 -c 0.95\'')
    subprocess.check_call(cmd, shell=True)

    os.remove(f"{outdir}/{prefix}.tmp.fasta")
    os.remove(f"{outdir}/{prefix}.merge.fasta")

if __name__=="__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("-p1", "--pe1", help="5' reads", required=True)
    parser.add_argument("-p2", "--pe2", help="3' reads", required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    parser.add_argument("-len", "--length", help="min assmbly length,default=150", type=int, default=150)
    args = parser.parse_args()
    run(args.pe1,args.pe2,args.prefix,args.outdir,args.length)