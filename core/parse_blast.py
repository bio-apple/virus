import os
import argparse
import subprocess

docker="virus:latest"
def run(blast_out_vsp,blast_out_nt_viruses,nt_virus_db_dir,outdir,accession):
    outfile = open(outdir + "/ref.list", "w")
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    db_dir = os.path.abspath(os.path.dirname(nt_virus_db_dir))
    db_name = nt_virus_db_dir.split("/")[-1]
    if accession is None:
        accession=[]
    else:
        for key in accession:
            outfile.write(f"{key}\n")
    num,query=len(accession),{}
    infile = open(blast_out_vsp, 'r')
    for line in infile:
        line = line.strip()
        array=line.split("\t")
        if not line.startswith("#") and array[0] in query:
            query[array[0]]=1
            if not array[1] in accession:
                accession.append(array[1])
                outfile.write(array[1] + "\n")
    infile.close()

    infile=open(blast_out_nt_viruses,'r')
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        if not line.startswith("#"):
            if not array[0] in query:
                num+=1
                if not array[1] in accession:
                    accession.append(array[1])
                    outfile.write(array[1]+"\n")
    infile.close()
    outfile.close()
    print(accession)
    # get species reference sequence from nt_viruses
    if num > 0:
        cmd = (f"docker run -v {db_dir}:/ref -v {outdir}:/outdir {docker} sh -c \'"
               f"export PATH=/opt/conda/envs/kraken2/bin/:$PATH && "
               f"export BLASTDB=/ref/ && "
               f"blastdbcmd -db /ref/{db_name} "
               f"-entry_batch /outdir/ref.list -out /outdir/ref.fasta\'")
        print(cmd)
        subprocess.run(cmd, shell=True)

    # get non-redundant fna and build reference index
    cmd = (f'docker run --rm -v {outdir}/:/outdir/ {docker} '
           f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
           f'cd-hit-est -i /outdir/ref.fasta -o /outdir/ref.non-redundant.fna -c 0.95 -M 0 && '
           f'bowtie2-build /outdir/ref.fasta /outdir/ref.fasta && '
           f'samtools faidx /outdir/ref.fasta && '
           f'bwa index -a bwtsw /outdir/ref.fasta\'')
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    print("Done")
    return num

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='blast2VSP')
    parser.add_argument("-v",'--blast_out_vsp',help='blast output file from vsp',required=True)
    parser.add_argument("-n","--blast_out_nt_viruses",help='blast output file from nt_viruses',required=True)
    parser.add_argument("-d","--nt_virus_dir",help='NT virus directory',required=True)
    parser.add_argument("-o","--outdir",help='output directory',default=os.getcwd())
    parser.add_argument("-a","--accession",help='accession',default=None)
    args = parser.parse_args()
    run(args.blast_out_vsp,args.blast_out_nt_viruses,args.nt_virus_dir,args.outdir,args.accession)


