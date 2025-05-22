import os,sys,re
import subprocess
import argparse

def run(bed,bam,outdir,prefix):
    bed=os.path.abspath(bed)
    bam=os.path.abspath(bam)
    outdir=os.path.abspath(outdir)
    if not os.path.exists(outdir):
        subprocess.check_call(f'mkdir -p {outdir}',shell=True)
    out=os.path.abspath(outdir)+"/"+prefix

    # trim primers with ivar (soft clipping)
    # https://andersen-lab.github.io/ivar/html/manualpage.html
    # -e    Include reads with no primers
    cmd = "ivar trim -e -i %s.bam -b %s -p %s.soft.clipped | tee %s.ivar.stdout && rm -rf %s.bam %s.bam.bai" % (out, bed, out, out, out, out)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    ## remove soft-clipped primers
    # https://jvarkit.readthedocs.io/en/latest/Biostar84452/
    # source activate && conda deactivate
    cmd = "samtools sort %s.soft.clipped.bam -o %s.soft.clipped.sort.bam" % (out, out)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    cmd = "java -jar jvarkit.jar biostar84452 --samoutputformat BAM %s.soft.clipped.sort.bam |samtools sort -n >%s.trimmed.bam" % (out, out)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    cmd = ("samtools fastq -1 %s_no_primer.R1.fq -2 %s_no_primer.R2.fq -s %s.singleton.fastq %s.trimmed.bam &>%s.bam2fastq.stdout"
           % (out, out, out, out, out))
    subprocess.check_call(cmd, shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Trim primers with ivar.")
    parser.add_argument("-m","--bam", help="sort bam file", required=True)
    parser.add_argument( "-b","--bed", help="reference fasta", required=True)
    parser.add_argument("-o","--outdir", help="output directory", required=True)
    parser.add_argument("-p","--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.bed, args.bam, args.outdir, args.prefix)
