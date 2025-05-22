import subprocess
import argparse

def scrub(R1,R2,database,prefix,outdir):
    out = outdir + "/" + prefix
    cmd = "cd %s && /software/seqtk/seqtk mergepe %s %s >%s.interleave.fq" % (outdir,R1,R2,prefix)
    subprocess.check_call(cmd, shell=True)
    cmd = ("cd %s && export PATH=/software/Python-v3.11.0/bin/:$PATH && "
           "/software/sra-human-scrubber/scripts/scrub.sh -s -p 24 -d %s -i %s.interleave.fq && "
           "/software/fastqtk/fastqtk deinterleave %s.interleave.fq.clean %s_no_host.R1.fq %s_no_host.R2.fq"%(outdir,database,prefix,prefix,prefix,prefix))
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call("rm -rf %s.interleave.fq %s.interleave.fq.clean"%(out,out),shell=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser("Human reads are removed with a modified version of the SRA Human Read Scrubber tool.")
    parser.add_argument("-1", "--R1", help="sort bam file", required=True)
    parser.add_argument("-2", "--R2", help="reference fasta", required=True)
    parser.add_argument("-d","--database",help="human k-mer database used by the NCBI SRA Human Read Removal Tool",required=True)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p","--prefix",help="prefix of output",required=True)
    args=parser.parse_args()
    scrub(args.R1,args.R2,args.database,args.prefix,args.outdir)