import os,sys,re
import subprocess
import argparse
import matplotlib.pyplot as plt
import numpy as np
import statistics
import pandas as pd

def run(bam,outdir,prefix):

    bam=os.path.abspath(bam)
    out=outdir+"/"+prefix
    result = subprocess.run(f"samtools view -H {bam}", shell=True, capture_output=True, text=True, check=True)
    if not "SO:coordinate" in result.stdout:
        sort_bam = os.path.splitext(bam)[0] + "sort.bam"
        subprocess.run(f"samtools sort -o {sort_bam} {bam}", shell=True, check=True)
    else:
        sort_bam = bam
    if not os.path.exists(sort_bam+".bai"):
        subprocess.check_call(f'samtools index {sort_bam}', shell=True)

    subprocess.check_call(f'samtools depth -J -d 8000 -Q 0 -q 20 -aa {sort_bam} >{out}.site.depth.txt',shell=True)
    infile = open("%s.depth.txt" % out, "r")
    cov = [0, 0]
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        if int(array[2]) == 0:
            cov[0] += 1
        if int(array[2]) < 10:
            cov[1] += 1
    infile.close()
    ##############plot coverage###########################
    df = pd.read_csv("%s.depth.txt" % (out),
                     sep='\t',
                     # engine='python',
                     names=["ref", "pos", "depth"]
                     )

    median_depth = statistics.median(df["depth"])
    plt.figure(figsize=(10, 4))
    plt.axhline(median_depth, linestyle='--', color='red', linewidth=1, label="median: %.0f" % median_depth)
    plt.axhline(10, linestyle='--', color='grey', linewidth=1, label="<10X(%s bp)" % (cov[1]))

    max = np.max(10000)
    maxlog10 = np.ceil(np.log10(max))
    plt.ylim(top=10 ** maxlog10)

    plt.title("Sample: %s\n" % (prefix), fontsize=10, wrap=True)
    plt.xlabel("Position along genome [bp]")
    plt.ylabel("Coverage depth")
    plt.yscale("log")
    plt.margins(x=0.01)
    plt.legend()
    plt.ylim(bottom=1)
    plt.yscale("log")
    plt.plot(df["pos"], df["depth"])
    plt.savefig("%s.coverage.png" % (out), dpi=300)
    print("\nPre-process Done.\n")
