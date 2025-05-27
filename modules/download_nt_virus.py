import os,sys,re
import subprocess
import argparse
import json
from datetime import date


def run(outdir):
    subprocess.check_call("wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses-nucl-metadata.json",shell=True)
    subprocess.check_call(f'mkdir -p {outdir}',shell=True)

    with open("./nt_viruses-nucl-metadata.json", "r") as f:
        data = json.load(f)
        verison=data['last-updated'].split("T")[0]
        files=data["files"]
        outfile = open(f'{outdir}/{verison}.sh',"w")
        for key in files:
            outfile.write(f"nohup wget {key}&\n")
        outfile.close()
    #subprocess.check_call(f'sh {outdir}/{verison}.sh',shell=True)

if __name__=="__main__":
    if len(sys.argv)!=2:
        print(f"\nusage:python3 {sys.argv[0]} /outdir/\n")
    else:
        run(sys.argv[1])