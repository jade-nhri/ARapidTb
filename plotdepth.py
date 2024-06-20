#!/usr/bin/env python3
import pandas as pd
import os,sys
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import Counter

indir=sys.argv[1]
indir=os.path.abspath(indir)

inbam=os.path.join(indir,'calls_to_draft.bam')
bamfile=pysam.AlignmentFile(inbam,'rb')



lengths = bamfile.lengths
#print (lengths)
nchromosome = len(lengths)
#print (nchromosome)

names=[]
for i in range(0,nchromosome):
    names.append(bamfile.get_reference_name(i))
names=sorted(names)

k=0
f=plt.figure(figsize=(9, 6))
for m in range (0,5):
    for n in range (0,3):
        xpos=[]
        ycov=[]
        #print (names[k])
        for pileupcolumn in bamfile.pileup(names[k],0,):
            ycov.append(pileupcolumn.n)
        for i in range(0,len(ycov)):
            xpos.append(i)
     
        #print (len(xpos))
        #print (len(ycov))

        p=plt.subplot2grid((5,3),(m,n))
        plt.ylim(0,200)
        plt.xlim(0,lengths[k])
        p.axes.xaxis.set_ticks([])
        plt.plot(xpos,ycov)
        plt.title(names[k],y=0.5)
        k+=1

    

f.savefig(os.path.join(indir,'depth.pdf'))

