import os,sys,re
import argparse
import subprocess


docker="virus:latest"
def run(blast_out,VSP_accession,VSP_ssname,nt_virus_dir,outdir):

    outdir=os.path.abspath(outdir)
    os.makedirs(outdir,exist_ok=True)

    nt_virus_dir=os.path.abspath(nt_virus_dir)

    #parse blast ssname
    infile = open(blast_out,"r")
    blast_ssname=[]
    for line in infile:
        line = line.strip()
        if not line.startswith("#"):
            if line.split("\t")[9] not in blast_ssname:
                blast_ssname.append(line.split("\t")[9])
    infile.close()
    print(blast_ssname)
    #parse VSP accession
    accession=[]
    infile = open(VSP_accession,"r")
    for line in infile:
        line=line.strip()
        accession.append(line)
    infile.close()
    print(accession)
    #output positive species accession
    outfile=open(outdir+"/ref.list","w")
    infile = open(VSP_ssname,"r")
    for line in infile:
        line=line.strip()
        array=line.split(",")
        if array[0].split(".")[0] in accession and array[1] in blast_ssname:
            outfile.write(array[0]+"\n")
    infile.close()
    outfile.close()
    #get species reference sequence from nt_viruses
    cmd=(f"docker run -v {nt_virus_dir}:/ref -v {outdir}:/outdir {docker} sh -c \'"
         f"export PATH=/opt/conda/envs/kraken2/bin/:$PATH && "
         f"export BLASTDB=/ref/ && "
         f"blastdbcmd -db /ref/nt_viruses "
         f"-entry_batch /outdir/ref.list -out /outdir/ref.fasta\'")
    print(cmd)
    subprocess.run(cmd,shell=True)

    #build reference index
    cmd=(f'docker run --rm -v {outdir}/:/outdir/ {docker} '
         f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
         f'bowtie2-build /outdir/ref.fasta /outdir/ref.fasta && '
        f'samtools faidx /outdir/ref.fasta && '
        f'bwa index -a bwtsw /outdir/ref.fasta\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    print("Done")

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='blast2VSP')
    parser.add_argument("-b",'--blast_out',help='blast output file',required=True)
    parser.add_argument("-a","--accession",help='VSP accession file',required=True)
    parser.add_argument("-s","--ssname",help='VSP ssname file',required=True)
    parser.add_argument("-n","--nt_virus_dir",help='NT virus directory',required=True)
    parser.add_argument("-o","--outdir",help='output directory',default=os.getcwd())
    args = parser.parse_args()
    run(args.blast_out,args.accession,args.ssname,args.nt_virus_dir,args.outdir)


