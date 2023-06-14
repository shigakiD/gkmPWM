import os,sys

def main():
    name = sys.argv[1]
    fidseq = open(name+'_svseq.fa', 'w')
    fidalpha = open(name+'_svalpha.out', 'w')
    fidmodel = open(name+'.model.txt','r')
    model = fidmodel.readlines()
    a = 1
    l = len(model)
    ind = 0
    for i in range(l):
        if ind == 1:
            line = model[i]
            x = line.split(' ')
            fidseq.write('>seq'+str(a)+'\n'+x[1].strip(' '))
            fidalpha.write('>seq'+str(a)+'\t'+x[0]+'\n')
            a = a+1
        if model[i].strip() == 'SV':
            ind = 1;
    fidseq.close()
    fidalpha.close()
    fidmodel.close()
main()
