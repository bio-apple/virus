import os
import argparse
import core
from concurrent.futures import ThreadPoolExecutor, as_completed


# 获取当前脚本的绝对路径
script_path = os.path.abspath(__file__)

# 获取脚本所在目录
script_dir = os.path.dirname(script_path)

parser=argparse.ArgumentParser("Virus NGS pipeline\nEamile:fanyucai3@gmail.com")
parser.add_argument("-p1","--pe1",help="R1 fastq",required=True, nargs='+')
parser.add_argument("-p2","--pe2",help="R2 fastq",default=None,nargs='+')
parser.add_argument("-p","--prefix",help="prefix of output",required=True, nargs='+')
parser.add_argument("-b","--bed",help="bed file",default=None)
parser.add_argument('-k','--kraken2',help='kraken2 reference index',required=True)
#parser.add_argument('-r','--ref',help="directory reference bowtie2 index",required=True)
#parser.add_argument('-f','--fna',help="fasta reference",required=True)
parser.add_argument("-host","--host",help="directory host bowtie2 index",required=True)
parser.add_argument("-o","--outdir",help="diretory of output",required=True)
args=parser.parse_args()

args.outdir=os.path.abspath(args.outdir)
os.makedirs(args.outdir,exist_ok=True)


for r1,r2,prefix in zip(args.pe1,args.pe2,args.prefix):
    print("""
        # ------------------------
        # Step 1: fastp qc
        # ------------------------
    """)
    #core.fastp.run(r1,args.outdir+"/1.fastp",prefix,r2)

    print("""
        # ------------------------
        # Step 2: kraken2
        # ------------------------
    """
    )
    #core.kraken2.run(r1,args.kraken2,prefix,args.outdir+"/2.kraken2",r2)

    print("""
        # ------------------------
        # Step 3: bowtie2 host filter
        # ------------------------
    """)
    #core.filter_host.run(r1,args.outdir+"/3.filter_host",args.host,prefix,r2)

    print("""
        # ------------------------
        # Step 4: denovo genome assembly:megahit and metaspades
        # ------------------------
    """)
    read1,read2="",""
    if r2:
        read1=args.outdir+"/"+"3.filter_host/"+prefix+"_1.fastq"
        read2=args.outdir+"/"+"3.filter_host/"+prefix+"_2.fastq"
    else:
        read1 = args.outdir + "/" + "3.filter_host/" + prefix +".unaligned.fastq"
        read2=None
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(core.megahit.run, read1, prefix, args.outdir + "/4.assembly/", read2),
            executor.submit(core.metaspades.run, read1, prefix, args.outdir + "/4.assembly/", read2)
        ]
        for future in as_completed(futures):
            print(future.result())


    # ------------------------
    # Step 3: bowtie2 mapping reference
    # ------------------------


      # ------------------------
      # Step 4: trim primer
      # ------------------------

  # ------------------------
  # Step 5: variant calling and consensus sequence
  # ------------------------

  # ------------------------
  # Step 6: plot coverage
  # ------------------------


  # ------------------------
  # Step 7: run nextclade and pangolin
  # ------------------------
