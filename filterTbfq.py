#!/usr/bin/env python3
import os,sys
import subprocess
import io,gzip
import pandas as pd

infile=sys.argv[1]
reffile=sys.argv[2]
outdir=os.path.abspath(sys.argv[3])
print (outdir)

if '.fastq.gz' in infile:
    tmp=os.path.split(infile)[-1]
    outfile=tmp.replace('.fastq.gz','_filtered.fastq') 
if '.fastq' in infile and '.gz' not in infile: 
    tmp=os.path.split(infile)[-1]
    outfile=tmp.replace('.fastq','_filtered.fastq')

comm='minimap2 {0} {1} -t 16 -w5 > out.paf'.format(reffile,infile)
print (comm)
subprocess.getoutput(comm)

def checkreads(df,fadict):
    dfset=df[df['Tname']=='rrs']
    dfsetin=dfset[(dfset['Tstart']<=300) | (dfset['Tend']>=1900)]
    dfsetout=dfset[(dfset['Tstart']>300) & (dfset['Tend']<1900)]
    temp=dfsetin['Qname'].values.tolist()
    readout=temp
    readid=dfsetout['Qname'].values
    for i in readid:
        tempseq=fadict[i][1]
        if 'GCAGACTAGAGTACT'in tempseq or 'AGTACTCTAGTCTGC' in tempseq or 'GACGCGTCTAGAGAT' in tempseq or 'ATCTCTAGACGCGTC' in tempseq or 'GACTCGTGAGAGACT' in tempseq or 'AGTCTCTCACGAGTC' in tempseq:
            #print(i)
            readout.append(i)

    dfset=df[df['Tname']=='rpoB']
    dfsetin=dfset[dfset['Tstart']<1900]
    dfsetout=dfset[dfset['Tstart']>=1900]
    temp=dfsetin['Qname'].values.tolist()
    readout.extend(temp)
    readid=dfsetout['Qname'].values
    for i in readid:
        tempseq=fadict[i][1]
        if 'GTGCGTGTGTATGTG' in tempseq or 'CACATACACACGCAC' in tempseq or 'TGTCGTGCACGCTGC' in tempseq or 'GCAGCGTGCACGACA' in tempseq:
            print (i)
            readout.append(i)

    #print (len(readout))
    return (readout)
    
    
d=dict()
if '.gz' in infile:
    print ('Reading fq in gz...')
    f=io.TextIOWrapper(gzip.open(infile,'r'))

else:
    f=open(infile)


while True:
    h=f.readline()
    if not h: break
    h=h.replace('\n','')
    #print (h)
    rID=h.split()[0]
    rID=rID.replace('@','')
    #print (rID)
    seq=f.readline().replace('\n','')
    qh=f.readline().replace('\n','')
    qual=f.readline().replace('\n','')
    d[rID]=[]
    d[rID].append(h)
    d[rID].append(seq)
    d[rID].append(qual)
f.close()

df=pd.read_table('out.paf',names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'])
print(df)
df=df[df['MQ1']>=20]
print(df)
df=df.drop_duplicates('Qname')
print(df)


df1=df[~df['Tname'].isin(['rrs','rpoB'])]
print (df1)
df2=df[df['Tname'].isin(['rrs','rpoB'])]
print (df2)

readID=df1['Qname'].values.tolist()
#print (readID)
readID2=checkreads(df2,d)
readID.extend(readID2)
print (len(readID))
#print (readID)

print(os.path.join(outdir,outfile))
fw=open(os.path.join(outdir,outfile),'w')
for i in readID:
    fw.write('@'+i+'\n')
    fw.write(d[i][1]+'\n')
    fw.write('+'+'\n')
    fw.write(d[i][2]+'\n')

fw.close()
