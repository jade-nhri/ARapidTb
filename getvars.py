#!/usr/bin/env python3
import os,sys
import pandas as pd
import pysam
from collections import Counter
from Bio.Seq import Seq
import json
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")



invcf=sys.argv[1] #vcf
infa=sys.argv[2] #TB_ref_amplicon
inbam=sys.argv[3]
varres='/opt/ARapidTb/Tbresdb.txt'

f=open(varres,'r')
dres=eval(f.read())

#variants not associated with resistance
dvarna={'rrs': {'C-187T': ['amikacin', 'capreomycin', 'kanamycin', 'linezolid', 'streptomycin'],'C517T': ['amikacin', 'capreomycin', 'kanamycin'], 'C492T': ['amikacin', 'capreomycin', 'kanamycin', 'streptomycin'], 'A906G': ['capreomycin']}, 'eis': {'V163I': ['amikacin', 'kanamycin'], 'C-12T': ['amikacin']}, 'embC': {'R738Q': ['ethambutol'], 'V981L': ['ethambutol'], 'C-565T': ['ethambutol'], 'T270I': ['ethambutol'], 'N394D': ['ethambutol'], 'C-1520A': ['ethambutol'], 'G-1743A': ['ethambutol'], 'C-900T': ['ethambutol'], 'A774S': ['ethambutol'], 'A-20C': ['ethambutol'], 'C-1193T': ['ethambutol'], 'V104M': ['ethambutol'], 'R567H': ['ethambutol'], 'L661I': ['ethambutol'], 'G-1419A': ['ethambutol'], 'C-589G': ['ethambutol'], 'G-270A': ['ethambutol'], 'C-4586G': ['ethambutol'], 'C-100T': ['ethambutol']}, 'ubiA': {'E149D': ['ethambutol'], 'V49I': ['ethambutol'], 'G268D': ['ethambutol'], 'A-3741C': ['ethambutol']}, 'embB': {'E378A': ['ethambutol'], 'N13S': ['ethambutol'], 'S1054P': ['ethambutol'], 'Q139H': ['ethambutol'], 'R213Q': ['ethambutol']}, 'katG': {'R463L': ['isoniazid'], 'C-85T': ['isoniazid'], 'C-354T': ['isoniazid'], 'V469L': ['isoniazid']}, 'ahpC': {'G-88A': ['isoniazid']}, 'inhA': {'C-40T': ['isoniazid'], 'C-522G': ['isoniazid']}, 'gyrA': {'E21Q': ['levofloxacin', 'moxifloxacin'], 'G668D': ['levofloxacin', 'moxifloxacin'], 'S95T': ['levofloxacin', 'moxifloxacin'], 'G247S': ['levofloxacin', 'moxifloxacin'], 'A384V': ['levofloxacin', 'moxifloxacin'], 'A463S': ['levofloxacin', 'moxifloxacin'], 'C-34T': ['levofloxacin', 'moxifloxacin'], 'Q613E': ['levofloxacin'], 'R252L': ['levofloxacin']}, 'gyrB': {'C-165T': ['levofloxacin', 'moxifloxacin'], 'M291I': ['levofloxacin', 'moxifloxacin'], 'P94L': ['levofloxacin', 'moxifloxacin'], 'A403S': ['levofloxacin'], 'V301L': ['levofloxacin', 'moxifloxacin']}, 'rpoB': {'C-61T': ['rifampicin'], 'A-261G': ['rifampicin'], 'V695L': ['rifampicin'], 'E250G': ['rifampicin'], 'I925V': ['rifampicin'], 'E639D': ['rifampicin'], 'S388L': ['rifampicin'], 'K944E': ['rifampicin'], 'N381H': ['rifampicin']}, 'rpsL': {'T-165C': ['streptomycin'], 'C-259T': ['streptomycin']}}

vard=dict()

df=pd.read_table(invcf,names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'],comment='#')
#print (df)
df.loc[(df['CHROM']=='gyrA') & (df['POS']==300) & (df['ALT']=='CC'),'ALT']='C'
#print (df)

##katG   934  .  GC  CA
##ethA   669  .  GC   G
##ethA  138 .   CA  A
##pncA  900 .   T   TG
##pncA  757 .   CXXXX   C
##pncA  842 .   ACCAGGGTG   A
##fabG1_inhA    709 .   G   A   #20240523
sindexs=df[((df['CHROM']=='ethA') & (df['POS']==669) &(df['REF']=='GC') &(df['ALT']=='G'))|((df['CHROM']=='katG') & (df['POS']==934) &(df['REF']=='GC') &(df['ALT']=='CA'))|((df['CHROM']=='ethA') &(df['POS']==138)&(df['REF']=='CA') &(df['ALT']=='C'))|((df['CHROM']=='pncA') & (df['POS']==900) &(df['REF']=='T') &(df['ALT']=='TG'))|((df['CHROM']=='pncA') & (df['POS']==757) &(df['ALT']=='C')) | ((df['CHROM']=='pncA') & (df['POS']==842) &(df['REF']=='ACCAGGGTG') &(df['ALT']=='A'))].index

#print (sindexs)
if sum((df['CHROM']=='ethA') & (df['POS']==669) &(df['REF']=='GC') &(df['ALT']=='G'))==1:
    vard['ethA']=['597delC']
if sum((df['CHROM']=='katG') & (df['POS']==934) &(df['REF']=='GC') &(df['ALT']=='CA'))==1:
    vard['katG']=['S315T']
if sum((df['CHROM']=='ethA') & (df['POS']==138) &(df['REF']=='CA') &(df['ALT']=='C'))==1:
    if 'ethA' in vard.keys():
        vard['ethA'].append('65delA')
    else:
        vard['ethA']=['65delA']
if sum((df['CHROM']=='pncA') & (df['POS']==900) &(df['REF']=='T') &(df['ALT']=='TG'))==1:
    vard['pncA']=['517dupG']
if sum((df['CHROM']=='pncA') & (df['POS']==842) &(df['REF']=='ACCAGGGTG') &(df['ALT']=='A'))==1:
    if 'pncA' in vard.keys():
        vard['pncA'].append('458delCCAGGGTG')
    else:
        vard['pncA']=['458delCCAGGGTG']
if sum((df['CHROM']=='pncA') & (df['POS']==757) &(df['ALT']=='C'))==1:
    if 'pncA' in vard.keys():
        vard['pncA'].append('373del')
    else:
        vard['pncA']=['373del']

if sum((df['CHROM']=='fabG1_inhA') & (df['POS']==709) &(df['REF']=='G') &(df['ALT']=='A'))==1:
    vard['inhA']=['G-154A']



for i in vard.keys():
    if i in dres.keys():
        for j in vard[i]:
            if j in dres[i].keys():
                if j=='597delC':
                    print ('***{0}\t{1}\t{2}'.format(i,669,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))
                if j=='65delA':
                    print ('***{0}\t{1}\t{2}'.format(i,138,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))
                if j=='S315T':
                    print ('***{0}\t{1}\t{2}'.format(i,934,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))
                if j=='517dupG':
                    print ('***{0}\t{1}\t{2}'.format(i,900,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))
                if j=='458delCCAGGGTG':
                    print ('***{0}\t{1}\t{2}'.format(i,842,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))
                if j=='373del':
                    print ('***{0}\t{1}\t{2}'.format(i,757,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))
                if j=='G-154A':
                    print ('***{0}\t{1}\t{2}'.format(i,709,j))
                    print ('***Resistant durgs: {0}'.format(dres[i][j]))



df=df.drop(sindexs)
#print (df)

bamfile=pysam.AlignmentFile(inbam,'rb')

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def getvarsite(gene,pos):
    aasite=0
    cdspos=0
    posout=0
    #ahpC: 337:994
    if gene=='ahpC':
        if pos<=942 and pos>=337:
            cdspos=pos-337+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+337
        #prompter
        if pos<337 and pos>=232:
            aasite=pos-337
            posout=pos
    #eis: 144:1352
    if gene=='eis':
        if pos<=1352 and pos>=144:
            cdspos=pos-144+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+144
        #promoter
        if pos<144:
            aasite=pos-144
            posout=pos
    #embB: -40:3101
    if gene=='embB':
        if pos<=3101:
            cdspos=pos+40
            aasite,posout=nt2aa(cdspos)
            posout=posout-39
    #embC: 196:
    if gene=='embC':
        cdspos=pos-196+1
        aasite,posout=nt2aa(cdspos)
        posout=posout+196
    #ethA: 75:1544
    if gene=='ethA':
        if pos<=1544 and pos>=75:
            cdspos=pos-75+1
            aasite,posout=nt2aa(cdspos)
            posout=pos+75

    #fabG1_inhA: 101:884; 863:1672
    if gene=='fabG1_inhA':
        if pos<=884 and pos>=101:
            gene='fabG1'
            cdspos=pos-101+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+101
        if pos<=1672 and pos>=863:
            gene='inhA'
            cdspos=pos-863+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+863
        #promoter
        if pos<101:
            gene='fabG1'
            aasite=pos-101
            posout=pos
    #gyrA: 17:
    if gene=='gyrA':
        if pos>=17:
            cdspos=pos-17+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+17
    #gyrB: -1002:1026
    if gene=='gyrB':
        if pos<=1026:
            cdspos=pos+1002
            aasite,posout=nt2aa(cdspos)
            posout=posout-1001

    #katG: -10:2213
    if gene=='katG':
        if pos<=2213:
            cdspos=pos+10
            aasite,posout=nt2aa(cdspos)
            posout=posout-9
    #pncA: 386:946
    if gene=='pncA':
        if pos<=946 and pos>=386:
            cdspos=pos-386+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+386
        #prompter
        if pos<386:
            aasite=pos-386
            posout=pos
  
    #rpoB: -58
    if gene=='rpoB':
        cdspos=pos+58
        aasite,posout=nt2aa(cdspos)
        posout=posout-57

    #rpsL: 108:482
    if gene=='rpsL':
        if pos<=482 and pos>=108:
            cdspos=pos-108+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+108

    #rrs: 334:1870
    #rrssite,posout
    if gene=='rrs':
        if pos<=1870 and pos>=334:
            aasite=pos-334+1
            posout=pos
        if pos<334:
            aasite=pos-334
            posout=pos

    #tlyA: 373:1179
    if gene=='tlyA':
        if pos<=1179 and pos>=373:
            cdspos=pos-373+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+373
    #ubiA: 246:1154
    if gene=='ubiA':
        if pos<=1154 and pos>=246:
            cdspos=pos-246+1
            aasite,posout=nt2aa(cdspos)
            posout=posout+246



    #print ('{0}\t{1}'.format(gene,cdspos))
    return(gene,aasite,posout)

def nt2aa(pos):
    pos=pos-1	#0-numbering
    aasite=(pos//3)+1
    posout=(aasite-1)*3	#0-numbering, codon start
    return(aasite,posout)


d=dict()
f=open(infa)
with f as fp:
    for name, seq in read_fasta(fp):
        name=name.replace('>','')
        d[name]=seq
f.close()

refd=dict()
altd=dict()
dout=dict()
for i in d.keys():
    dfset=df[df['CHROM']==i]
    targetposs=dfset['POS'].values
    alt=dfset['ALT'].values
    ref=dfset['REF'].values
    pos=0
    for j,k,l in zip(targetposs,alt,ref):
        #print ('{0}\t{1}\t{2}'.format(j,k,l))
        if i not in altd.keys():
            altd[i]=dict()
            refd[i]=dict()
        altd[i][j]=k
        refd[i][j]=l
    for pileupcolumn in bamfile.pileup(i,0,):
        pos+=1
        bases=[]
        if pileupcolumn.n>=20 and pileupcolumn.pos+1 in targetposs:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del:
                    bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
            if len(Counter(bases))>=1:
                if i not in dout.keys():
                    dout[i]=dict()
                dout[i][pileupcolumn.pos+1]=dict(Counter(bases))
#print ('======Gene with variant positions======')
#print (altd)
for i in dout.keys():
    for j in dout[i].keys():
        if 'A' not in dout[i][j].keys():
            dout[i][j]['A']=0
        if 'T' not in dout[i][j].keys():
            dout[i][j]['T']=0
        if 'C' not in dout[i][j].keys():
            dout[i][j]['C']=0
        if 'G' not in dout[i][j].keys():
            dout[i][j]['G']=0

#print ('======Read depth======')
#print (dout)

#vard=dict()
for i in dout.keys():
    #print(i)
    for j in dout[i].keys():
        #print (j)
        l,m,n=getvarsite(i,int(j))
        #print ('{0}\t{1}\t{2}'.format(l,m,n))
        if l!='rrs' and n!=0 and m>0:
            codonref=d[i][n-1:n+2]
            #print (codonref)
            #print (Seq(codonref).translate())
            aa1=Seq(codonref).translate()
            codonvar=codonref
            if (j-n)==0:
                codonvar=altd[i][j]+d[i][n:n+2]
            if (j-n)==1:
                codonvar=d[i][n-1]+altd[i][j]+d[i][n+1]
            if (j-n)==2:
                codonvar=d[i][n-1:n+1]+altd[i][j]
            #print (codonvar)
            aa2=Seq(codonvar).translate()
            if aa1!=aa2:
                if l not in vard.keys():
                    vard[l]=[]
                vard[l].append(str(aa1)+str(m)+str(aa2))
                if l in dres.keys():
                    if str(aa1)+str(m)+str(aa2) in dres[l].keys():
                        print ('***{0}\t{1}\t{2}'.format(i,j,str(aa1)+str(m)+str(aa2)))
                        print ('***Resistant drugs:{0} with variant depth {1} compared to reference depth {2}'.format(dres[l][str(aa1)+str(m)+str(aa2)],dout[i][j][altd[i][j]],dout[i][j][refd[i][j]]))
                if l in dvarna.keys():
                    if str(aa1)+str(m)+str(aa2) in dvarna[l].keys():
                        print ('###{0}\t{1}\t{2}'.format(i,j,str(aa1)+str(m)+str(aa2)))
                        print ('###Variants not associated with resistance: {0} with variant depth {1} compared to reference depth {2}'.format(str(aa1)+str(m)+str(aa2),dout[i][j][altd[i][j]],dout[i][j][refd[i][j]]))


        if m<0:
            if l not in vard.keys():
                vard[l]=[]
            vard[l].append(str(d[i][n-1])+str(m)+str(altd[i][j]))
            if l in dres.keys():
                if str(d[i][n-1])+str(m)+str(altd[i][j]) in dres[l].keys():
                    print ('***{0}\t{1}\t{2}'.format(i,j,str(d[i][n-1])+str(m)+str(altd[i][j])))
                    print ('***Resistant drugs:{0} with variant depth {1} compared to reference depth {2}'.format(dres[l][str(d[i][n-1])+str(m)+str(altd[i][j])],dout[i][j][altd[i][j]],dout[i][j][refd[i][j]]))
            if l in dvarna.keys():
                if str(d[i][n-1])+str(m)+str(altd[i][j]) in dvarna[l].keys():
                    print ('###{0}\t{1}\t{2}'.format(i,j,str(d[i][n-1])+str(m)+str(altd[i][j])))
                    print ('###Variants not associated with resistance: {0} with variant depth {1} compared to reference depth {2}'.format(str(d[i][n-1])+str(m)+str(altd[i][j]),dout[i][j][altd[i][j]],dout[i][j][refd[i][j]]))

           
        if l=='rrs' and m>0 and n!=0:
            if l not in vard.keys():
                vard[l]=[]
            vard[l].append(str(d[i][n-1])+str(m)+str(altd[i][j]))
            if l in dres.keys():
                if str(d[i][n-1])+str(m)+str(altd[i][j]) in dres[l].keys():
                    print ('***{0}\t{1}\t{2}'.format(i,j,str(d[i][n-1])+str(m)+str(altd[i][j])))
                    print ('***Resistant drugs:{0} with variant depth {1} compared to reference depth {2}'.format(dres[l][str(d[i][n-1])+str(m)+str(altd[i][j])],dout[i][j][altd[i][j]],dout[i][j][refd[i][j]]))

            if l in dvarna.keys():
                if str(d[i][n-1])+str(m)+str(altd[i][j]) in dvarna[l].keys():
                    print ('###{0}\t{1}\t{2}'.format(i,j,str(d[i][n-1])+str(m)+str(altd[i][j])))
                    print ('###Variants not associated with resistance: {0} with variant depth {1} compared to reference depth {2}'.format(str(d[i][n-1])+str(m)+str(altd[i][j]),dout[i][j][altd[i][j]],dout[i][j][refd[i][j]]))


print ('======Variant positions======')
print (vard)
