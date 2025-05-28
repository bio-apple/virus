import os
import subprocess
import argparse

docker="virus:latest"

def run(pe1,prefix,outdir,pe2=None):

    ############prepare input##########################
    pe1=os.path.abspath(pe1)
    os.makedirs(outdir,exist_ok=True)
    ##############################################################################
    cmd=f'docker run --rm -v {pe1}:/raw_data/{pe1.split("/")[-1]} -v {outdir}:/outdir/ '
    if pe2:
        cmd+=(f'-v {pe2}:/raw_data/{pe2.split("/")[-1]} {docker} '
              f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
              f'metaspades.py --only-assembler --memory 500 --threads 48 -1 /raw_data/{pe1.split("/")[-1]} -2 /raw_data/{pe2.split("/")[-1]} -o /outdir/spades_{prefix}/ && '
              f'quast.py --plots-format png --min-contig 500 --no-html --no-icarus --output-dir /outdir/spades_{prefix}/ /outdir/spades_{prefix}/scaffolds.fasta\'')
    else:
        cmd += (f'{docker} '
                f'sh -c \'export PATH=/opt/conda/bin:$PATH && '
                f'metaspades.py --only-assembler --memory 500 --threads 48 -s /raw_data/{pe1.split("/")[-1]} -o /outdir/spades_{prefix}/ && '
                f'quast.py --plots-format png --min-contig 500 --no-html --no-icarus --output-dir /outdir/spades_{prefix}/ /outdir/spades_{prefix}/scaffolds.fasta\'')

    print(cmd)
    subprocess.check_call(cmd,shell=True)

    infile=open(f"{outdir}/spades_{prefix}/scaffolds.fasta","r")
    outfile1=open(f"{outdir}/spades_{prefix}/scaffolds_500bp.fasta","w")
    outfile2=open(f"{outdir}/spades_{prefix}/scaffolds_1000bp.fasta","w")
    outfile3=open(f"{outdir}/spades_{prefix}/scaffolds_1500bp.fasta","w")
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
    return "Genome assembly metaspades done."

if __name__=="__main__":
    parser = argparse.ArgumentParser("assemble genome using metaSPAdes.")
    parser.add_argument("-p1", "--pe1", help="5' reads", required=True)
    parser.add_argument("-p2", "--pe2", help="3' reads", default=None)
    parser.add_argument("-o", "--outdir", help="output directory", required=True)
    parser.add_argument("-p", "--prefix", help="prefix of output", required=True)
    args = parser.parse_args()
    run(args.pe1,args.prefix,args.outdir,args.pe2)

