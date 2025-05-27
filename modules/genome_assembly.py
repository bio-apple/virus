import os,sys,re
import argparse
import subprocess


def run(pe1,pe2,outdir,prefix,min_length=1000):
    spades=""
    megahit=""