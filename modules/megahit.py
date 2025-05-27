import os,sys,re
import subprocess
import argparse

docker="meta:latest"
def run(pe1,prefix,outdir,pe2=None):
    array=pe1.split(",")
    in_dir=os.path.dirname(os.path.abspath(array[0]))
    a="/raw_data/"+os.path.abspath(array[0]).split("/")[-1]
    outdir =os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        if os.path.exists(outdir+"/megahit_%s/"%(prefix)):
            print("Please check the output directory already exists")
            exit(1)
    cmd = "docker run -v %s:/raw_data/ -v %s:/outdir/ %s sh -c \'export PATH=/opt/conda/envs/metawrap/bin:$PATH && megahit" % (in_dir,outdir, docker)
    for i in range(1,len(array)):
        if in_dir!=os.path.dirname(os.path.abspath(array[i])):
            print("All sample reads must be in the same directory.")
            exit(1)
        a+=","+"/raw_data/"+os.path.abspath(array[i]).split("/")[-1]
    if pe2 is not None:
        b=""
        array = pe2.split(",")
        for i in range(0,len(array)):
            if in_dir!=os.path.dirname(os.path.abspath(array[i])):
                print("All sample reads must be in the same directory.")
                exit(1)
            if i==0:
                b="/raw_data/"+os.path.abspath(array[0]).split("/")[-1]
            else:
                b+=","+"/raw_data/"+os.path.abspath(array[i]).split("/")[-1]
        cmd+=" -1 %s -2 %s "%(a,b)
    else:
        cmd+=" --read %s"%(a)
    cmd+=("-o /outdir/megahit_%s/ --out-prefix %s -m 0.95 --continue --min-contig-len 1000 -t 64 && "
          "quast.py --threads 36 --plots-format png --no-html --no-icarus "
          "--output-dir /outdir/megahit_%s/ /outdir/megahit_%s/%s.contigs.fa\'")%(prefix,prefix,prefix,prefix,prefix)
    print(cmd)
    subprocess.call(cmd,shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser("MEGAHIT is an ultra-fast and memory-efficient NGS assembler")
    parser.add_argument("-p1", "--pe1", help="comma-separated list of fasta/q paired-end #1 files", required=True)
    parser.add_argument("-p2", "--pe2", help="comma-separated list of fasta/q paired-end #2 files")
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.pe1,args.prefix,args.outdir,args.pe2)
