import sys
import time
import os,re
import subprocess

docker="fanyucai1/virus:latest"
def run(db,query,outdir,prefix):
    outdir=os.path.abspath(outdir)
    db = os.path.abspath(db)
    query=os.path.abspath(query)
    out=outdir+"/"+prefix
    model="a+"
    if not os.path.exists("%s.nextclade.tsv" % (out)):
        model="w"
    outfile = open("%s.nextclade.tsv" % (out), model)
    if model == "w":
        outfile.write("seqName\tclade\tqc.overallStatus\n")#https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/04-results-tsv.html
    for root, dirs, files in os.walk(db):
        for dir in dirs:
            db_dir=db+"/"+dir
            cmd = ("docker run -v %s:/database/ -v %s:/mnt %s "
                   "nextclade run --silent -C seqName,clade,qc.overallStatus -D /database/ "
                   "--output-tsv /mnt/%s.%s.nextclade.tsv /mnt/%s") %(db_dir,outdir,docker,prefix,dir,query.split("/")[-1])
            subprocess.call(cmd,shell=True)
            infile=open("%s.%s.nextclade.tsv"%(out,dir),"r")
            for line in infile:
                line=line.strip()
                array=line.split("\t")
                if len(array)==3 and not re.search(r'unassigned',line) and array[0]!="seqName" and re.search(r'\S',array[1]):
                    outfile.write("%s\n"%line)
            os.remove("%s.%s.nextclade.tsv"%(out,dir))
    outfile.close()