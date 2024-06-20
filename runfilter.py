#!/usr/bin/env python3
import os,sys
import subprocess

indir=sys.argv[1]
reffile=sys.argv[2]
outdir=os.path.abspath(sys.argv[3])

cwd=os.getcwd()
os.chdir(indir)

myfiles=[x for x in os.listdir() if '.fastq' in x and 'filter' not in x]
print(myfiles)
os.chdir(cwd)

for i in sorted(myfiles):
    comm='filterTbfq.py {0}/{1} {2} {3}'.format(indir,i,reffile,outdir)
    print(comm)
    subprocess.getoutput(comm)
