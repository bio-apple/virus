import os
import subprocess
import argparse

docker="virus:latest"
parser=argparse.ArgumentParser("assemble genome using metaSPAdes.")
parser.add_argument("-p1","--pe1",help="5' reads",required=True)
parser.add_argument("-p2","--pe2",help="3' reads",required=True)
parser.add_argument("-o","--outdir",help="output directory",required=True)
parser.add_argument("-p","--prefix",help="prefix of output",required=True)
args=parser.parse_args()

############prepare input##########################
args.pe1=os.path.abspath(args.pe1)
args.pe2 = os.path.abspath(args.pe2)
in_dir=os.path.dirname(args.pe1)
if in_dir!=os.path.dirname(args.pe2):
    print("read1 and reads2 must be in the same directory.")
    exit()
a=args.pe1.split("/")[-1]
b=args.pe2.split("/")[-1]
args.outdir=os.path.abspath(args.outdir)
if not os.path.exists(args.outdir):
    subprocess.check_call('mkdir -p %s'%(args.outdir),shell=True)
##############################################################################
cmd=(f"docker run -v {in_dir}:/raw_data/ -v {args.outdir}:/outdir/ {docker} "
     f"sh -c \'export PATH=/opt/conda/bin:$PATH && "
     f"metaspades.py --only-assembler --memory 500 "
     f"--threads 48 -1 /raw_data/{a} -2 /raw_data/{b} -o /outdir/spades_{args.prefix}/")
cmd+=(f" && quast.py --plots-format png --min-contig 500 --no-html --no-icarus "
      f"--output-dir /outdir/spades_{args.prefix}/ /outdir/spades_{args.prefix}/scaffolds.fasta\'")
print(cmd)
subprocess.check_call(cmd,shell=True)

infile=open(f"{args.outdir}/spades_{args.prefix}/scaffolds.fasta","r")
outfile1=open(f"{args.outdir}/spades_{args.prefix}/scaffolds_500bp.fasta","w")
outfile2=open(f"{args.outdir}/spades_{args.prefix}/scaffolds_1000bp.fasta","w")
outfile3=open(f"{args.outdir}/spades_{args.prefix}/scaffolds_1500bp.fasta","w")
fa,id={},""
for line in infile:
    line=line.strip()
    if line.startswith(">"):
        id=line
        fa[id]=""
    else:
        fa[id]+=line
infile.close()
for key in fa:
    if int(key.split("_")[3])>=500:
        outfile1.write("%s\n%s\n"%(key,fa[key]))
    if int(key.split("_")[3]) >= 1000:
        outfile2.write("%s\n%s\n" % (key, fa[key]))
    if int(key.split("_")[3]) >= 1500:
        outfile3.write("%s\n%s\n" % (key, fa[key]))
outfile1.close()
outfile2.close()
outfile3.close()

