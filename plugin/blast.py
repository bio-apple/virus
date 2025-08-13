import os,sys,re
import subprocess,argparse

blast=""

def acc2name(accessionlist):
    infile = open(accessionlist,'r')
    for line in infile:
        line = line.strip()
        cmd="blastdbcmd -db nt -dbpath /path/to/blast_dbs -entry NC_034556 -outfmt \"%t\""
        subprocess.check_call(cmd,shell=True)
    infile.close()
