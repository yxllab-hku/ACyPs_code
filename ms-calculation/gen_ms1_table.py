import argparse
import pandas as pd
import os


working_path = str(os.getcwd())

etoga = argparse.ArgumentParser(description='example:python gen_ms1_table.py -seq GMYRPYIAKYVEEQTLQNSTNLVYDDITQISFINKEKNVKKINLGPDTTIVTETIENADPDEYFLDPDTTIHTFTVENNDTDEHYCGPETTRITKTLENIDPDEYYS -mod H2O -mod_num 6 -charge 20 -mn 1 -o test')
etoga.add_argument('-seq', type=str, nargs=1, required=True, help='input sequence')
etoga.add_argument('-mod',type=str,nargs=1,required=True, help='modification type,such as -H,-H2O> -mod_num <numbers of modification, range from 0 to very large number, you could estimate the number based on the residues involved in the modification')
etoga.add_argument('-mod_num',type=str,nargs=1, required=True,help='numbers of modification, range from 0 to very large number, you could estimate the number based on the residues involved in the modification')
etoga.add_argument('-charge',type=str,nargs=1, required=True, help='numbers of charge,range from 0 to very large number')
etoga.add_argument('-mn',type=str,nargs=1,required=True, help='the number of modification at the same time, such as dropping 2 H')
etoga.add_argument('-o', type=str, nargs=1, required=True, help='output filename, without suffix')

args = etoga.parse_args()


monomass={'A':71.037113805,'R':156.10111105,'N':114.04292747,'D':115.026943065,'C':103.009184505,'E':129.042593135,'Q':128.05857754,'G':57.021463735,'H':137.058911875,'I':113.084064015,'L':113.084064015,'K':128.09496305,'M':131.040484645,'F':147.068413945,'P':97.052763875,'S':87.032028435,'T':101.047678505,'W':186.07931298,'Y':163.063328575,'V':99.068413945}


def gen_multiseq(seq):
    seqlist=[]
    for i in range(len(seq)):
        for j in range(i,len(seq)):
            seqlist.append(seq[(len(seq)-j-1):(len(seq)-i+1)])
    return seqlist


def isomass(multiseq,monomass):
    im=[]
    for i in multiseq:
        m=0
        for j in i:
            m+=monomass[j]
        im.append(m+2*1.007825035+15.99491463)
    return im
       

def chargecolumn(im,cha):
    chargemass=[]
    for i in im:
        chargemass.append(float((i+cha*1.00727646688)/cha))
    return chargemass

def modify(im,mdn,mdv,cha):
    modmass=[]
    for i in im:
        modmass.append(float((i+mdn*mdv+cha*1.00727646688)/cha))
    return modmass

data={}
multiseq=gen_multiseq(str(args.seq[0]))
im=isomass(multiseq,monomass)
data['seq']=multiseq
data['isomass']=im
if str(args.mod[0])=='H2O':
   mdv=-(2*1.007825035+15.99491463)
if str(args.mod[0])=='H':
   mdv=-2*1.007825035
for i in range(1,int(args.charge[0])+1):
    data['+%s'%str(i)]=chargecolumn(im,i)

for j in range(1,int(args.mod_num[0])+1):
    for k in range(1,int(args.charge[0])+1):
        data['-%s%s+%s'%(str(j*int(args.mn[0])),str(args.mod[0]),k)]=modify(im,j,mdv,k)
           

df=pd.DataFrame(data)
df.to_csv('%s.csv'%str(args.o[0]),index=False)
 
