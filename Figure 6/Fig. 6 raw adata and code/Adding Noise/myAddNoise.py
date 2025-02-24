import sys
import numpy

CSV='C:/Users/tuo47278/Desktop/noise_maker_WD/PD9478_WildMutCount_Gene.csv'
RepNum=1
VAFcut=0
TotRead=30

In=open(CSV,'r').readlines()[1:]
Rep=1
while Rep<=RepNum:
    out=['\t'.join(['mutation_id','ref_counts','var_counts','normal_cn','minor_cn','major_cn'])+'\n']
    for i in In:
        i=i.split(',')
        ID=i[0]
        MutCell=int(i[2])
        WildCell=int(i[1])		
        VAF=1.0*MutCell/((MutCell+WildCell)*2)
        if VAF>VAFcut:
             SimTotRead=numpy.random.poisson(TotRead)
             Mut=numpy.random.binomial(SimTotRead,VAF)
             Ref=SimTotRead-Mut
            # print (ID,MutCell,WildCell,VAF,SimTotRead,Mut)
             out.append('\t'.join([ID,str(Ref),str(Mut),'2','0','2'])+'\n')
    Out=CSV[:-4]+'_D'+str(TotRead)+'_'+str(Rep)+'PyCloneIn.tsv'
    OutF=open(Out,'w')
    OutF.write(''.join(out))
    OutF.close()
    Rep+=1    
