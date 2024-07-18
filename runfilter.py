#!/usr/bin/env python3
import os,sys
import subprocess
import argparse

outpath='output'
refpath='/opt/ARapidTb/Tbrefs.fasta'
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='the path to fastq folder')
parser.add_argument('-r', help='the file path of target genes')
parser.add_argument('-o', help='an output folder (dfault output)')

args = parser.parse_args()

argv=sys.argv
if '-i' in argv:
    inpath=os.path.abspath(argv[argv.index('-i')+1])
if '-r' in argv:
    refpath=os.path.abspath(argv[argv.index('-r')+1])
if '-o' in argv:
    outpath=os.path.abspath(argv[argv.index('-o')+1])

indir=inpath
reffile=refpath
outdir=outpath
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
cwd=os.getcwd()
os.chdir(indir)

myfiles=[x for x in os.listdir() if '.fastq' in x and 'filter' not in x]
print(myfiles)
os.chdir(cwd)

for i in sorted(myfiles):
    comm='filterTbfq.py {0}/{1} {2} {3}'.format(indir,i,reffile,outdir)
    print(comm)
    subprocess.getoutput(comm)
