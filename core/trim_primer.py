import os,sys,re
import subprocess
import argparse

docker="virus:latest"
def run(bed,bam,outdir,prefix):
    bed=os.path.abspath(bed)
    bam=os.path.abspath(bam)
    outdir=os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)

    cmd=(f'docker run --rm '
         f'-v {bed}:/raw_data/{bed.split("/")[-1]} '
         f'-v {bam}:/raw_data/{bam.split("/")[-1]} '
         f'-v {outdir}:/outdir/ {docker} sh -c \''
         f'export PATH=/opt/conda/bin:$PATH && export JAVA_HOME=/opt/conda/ && ')

    # trim primers with ivar (soft clipping)
    # https://andersen-lab.github.io/ivar/html/manualpage.html
    # -e    Include reads with no primers
    ivar =cmd+ (f'ivar trim -e -i /raw_data/{bam.split("/")[-1]} -b /raw_data/{bed.split("/")[-1]} -p /outdir/{prefix}.soft.clipped '
            f'| tee /outdir/{prefix}.ivar.stdout && rm -rf /outdir/{prefix}.bam /outdir/{prefix}.bam.bai\'')
    print(ivar)
    subprocess.check_call(ivar,shell=True)

    ## remove soft-clipped primers
    # https://jvarkit.readthedocs.io/en/latest/Biostar84452/
    # source activate && conda deactivate
    out_bam=cmd+ f"samtools sort /outdir/{prefix}.soft.clipped.bam -o /outdir/{prefix}.soft.clipped.sort.bam && jvarkit biostar84452 --samoutputformat BAM /outdir/{prefix}.soft.clipped.sort.bam |samtools sort -n >/outdir/{prefix}.trimmed.bam\'"
    print(out_bam)
    subprocess.check_call(out_bam,shell=True)

    out_fastq=cmd+f"samtools fastq -1 /outdir/{prefix}_no_primer.R1.fq -2 /outdir/{prefix}_no_primer.R2.fq -s /outdir/{prefix}.singleton.fastq /outdir/{prefix}.trimmed.bam &>/outdir/{prefix}.bam2fastq.stdout\'"
    print(cmd)
    subprocess.check_call(out_fastq, shell=True)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Trim primers with ivar.")
    parser.add_argument("-m","--bam", help="sort bam file", required=True)
    parser.add_argument( "-b","--bed", help="reference fasta", required=True)
    parser.add_argument("-o","--outdir", help="output directory", required=True)
    parser.add_argument("-p","--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.bed, args.bam, args.outdir, args.prefix)
