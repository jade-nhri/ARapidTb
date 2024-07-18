#!/usr/bin/env python3
import os,sys
import subprocess
import pandas as pd
import argparse

outpath='output'
refpath='/opt/ARapidTb/Tbrefs.fasta'
parser = argparse.ArgumentParser()
parser.add_argument('-i', help='the path to the filtered fastq folder')
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
ref=refpath
outdir=outpath

if not os.path.exists(outdir):
    os.mkdir(outdir)

cwd=os.getcwd()
print (cwd)

os.chdir(indir)
myfiles=[x for x in os.listdir() if '_filtered' in x]
#print (myfiles)
os.chdir(cwd)

for i in sorted(myfiles):
    name=i.replace('_filtered.fastq','')
    print (name)
    comm='medaka_consensus -i {0}/{1} -d {2} -o {3}/medaka_{4}'.format(indir,i,ref,outdir,name)
    #print (comm)
    subprocess.getoutput(comm)

    comm='samtools idxstats {0}/medaka_{1}/calls_to_draft.bam > {0}/medaka_{1}/readcounts.txt'.format(outdir,name)
    #print (comm)
    subprocess.getoutput(comm)

    rc=pd.read_table('{0}/medaka_{1}/readcounts.txt'.format(outdir,name),names=['gene','length','count','comm'])
    rcset=rc[rc['gene']!='*']
    rcset=rcset[rcset['count']>=40]
    if len(rcset)>=3:
        print (rcset['gene'].values)
        comm='medaka variant {2} {0}/medaka_{1}/consensus_probs.hdf {0}/variant_{1}.vcf'.format(outdir,name,ref)
        #print (comm)
        stdout=subprocess.getoutput(comm)
        #print (stdout)

        comm='getvars.py {0}/variant_{1}.vcf {2} {0}/medaka_{1}/calls_to_draft.bam'.format(outdir,name,ref)
        stdout=subprocess.getoutput(comm)
        print (stdout)
    comm='plotdepth.py {0}/medaka_{1}'.format(outdir,name)
    subprocess.getoutput(comm)
    print ('\n')
