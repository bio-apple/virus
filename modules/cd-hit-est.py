# Eamil:yucai.fan@illumina.com
import os
import subprocess
import argparse
import time

parser=argparse.ArgumentParser("")
parser.add_argument("-f","--fna",help="fasta sequence",required=True)
parser.add_argument("-c","--identify",help="sequence identity threshold, default: 0.998",default=0.998,type=float,required=True)
parser.add_argument("-p","--prefix",help="prefix of output",default=time.strftime("%Y-%m-%d"))
parser.add_argument("-o","--outdir",help="output directory",default=os.getcwd(),required=True)
args=parser.parse_args()

start=time.time()
args.outdir=os.path.abspath(args.outdir)
args.fna=os.path.abspath(args.fna)
if not os.path.exists(args.outdir):
    subprocess.check_call('mkdir -p %s'%args.outdir,shell=True)
# Lythgoe K A, Hall M, Ferretti L, et al. SARS-CoV-2 within-host diversity and transmission[J]. Science, 2021, 372(6539): eabg0821.
# 0.995 sequence identity
# Dezordi F Z, Resende P C, Naveca F G, et al. Unusual SARS-CoV-2 intrahost diversity reveals lineage superinfection[J]. Microbial Genomics, 2022, 8(3).
# 99.8% sequence identity
out=args.outdir+"/"+args.prefix
cmd="cd-hit-est -i %s -o %s.fna -c %s"%(args.fna,out,args.identify)
subprocess.check_call(cmd,shell=True)

end=time.time()
print("Elapse time is %g seconds" % (end - start))