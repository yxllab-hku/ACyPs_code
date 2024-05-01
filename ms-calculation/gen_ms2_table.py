import argparse
import pandas as pd
import os
import sys
import numpy as np

working_path = str(os.getcwd())

etoga = argparse.ArgumentParser(description='example:python gen_ms2_table.py -seq GMYRPYIAKYVEEQTLQNSTNLVYDDITQISFINKEKNVKKINLGPDTTIVTETIENADPDEYFLDPDTTIHTFTVENNDTDEHYCGPETTRITKTLENIDPDEYYS -mod H2O -mod_num 2 -charge 10 -mn 1 -o test')
etoga.add_argument('-seq', type=str, nargs=1, required=True, help='input sequence')
etoga.add_argument('-mod',type=str,nargs=1, help='modification type,such as -H,-H2O> -mod_num <numbers of modification, range from 0 to very large number, you could estimate the number based on the residues involved in the modification')
etoga.add_argument('-mod_num',type=str,nargs=1,help='numbers of modification, range from 0 to very large number, you could estimate the number based on the residues involved in the modification')
etoga.add_argument('-charge',type=str,nargs=1, required=True, help='numbers of charge,range from 0 to very large number')
etoga.add_argument('-mn',type=str,nargs=1, help='the number of modification at the same time, such as dropping 2 H')
etoga.add_argument('-o', type=str, nargs=1, required=True, help='output filename, without suffix')

args = etoga.parse_args()


if args.mod is not None and (args.mod_num is None or args.mn is None):
   print('please provide the value of mod_num and mn!!!')
   sys.exit()

monomass={'A':71.037113805,'R':156.10111105,'N':114.04292747,'D':115.026943065,'C':103.009184505,'E':129.042593135,'Q':128.05857754,'G':57.021463735,'H':137.058911875,'I':113.084064015,'L':113.084064015,'K':128.09496305,'M':131.040484645,'F':147.068413945,'P':97.052763875,'S':87.032028435,'T':101.047678505,'W':186.07931298,'Y':163.063328575,'V':99.068413945}

def gen_multiseq(seq,sig,mod,mod_num,mn):
    mindex=[]
    seqlist=[]
    for i in range(len(seq)):
        if mod==None:
           mindex.append('%s%s'%(sig,str(i+1)))
        else:
           mindex.append('%s%s-%s%s'%(sig,str(i+1),str(mod_num*mn),str(mod)))
        if sig=='b':
           seqlist.append(seq[:(i+1)])
        else:
           seqlist.append(seq[(len(seq)-i-1):])
    return mindex,seqlist


def isomass(multiseq,monomass,sig,mod,mod_num,mn):
    im=[]
    if mod is None:
       mdv=0
    if mod=='H':
       mdv=-mn*mod_num*1.007825035
    if mod=='H2O':
       mdv=-mn*mod_num*(2*1.007825035+15.99491463)
    if mod=='CH3':
       mdv=mn*mod_num*14.01565007
    for i in multiseq:
        m=0
        for j in i:
            m+=monomass[j]
        if sig=='y':
           im.append(m+2*1.007825035+15.99491463+mdv)
        else:
           im.append(m+mdv)
    return im


def chargecolumn(im,cha):
    chargemass=[]
    for i in im:
        chargemass.append(float((i+cha*1.00727646688)/cha))
    return chargemass

ms2=['y','b']
writer=pd.ExcelWriter('%s.xlsx'%str(args.o[0]))
for ms in ms2:
    data={}
    data['index'],data['seq']=gen_multiseq(str(args.seq[0]),ms,None,None,None)
    im=isomass(data['seq'],monomass,ms,None,None,None)
    for c in range(1,int(args.charge[0])+1):
        data['ion+%s'%str(c)]=chargecolumn(im,c)
    if args.mod is not None:
       for num in range(1,int(args.mod_num[0])+1):
           multindex,multiseq=gen_multiseq(str(args.seq[0]),ms,str(args.mod[0]),num,int(args.mn[0]))
           data['index']+=multindex
           data['seq']+=multiseq
           im=isomass(multiseq,monomass,ms,str(args.mod[0]),num,int(args.mn[0]))
           for c in range(1,int(args.charge[0])+1):
               data['ion+%s'%str(c)]+=chargecolumn(im,c)
    df=pd.DataFrame(data)
    df.to_excel(writer,sheet_name='%s_ion'%ms,index=False)
writer.save()


