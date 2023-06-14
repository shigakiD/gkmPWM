import logomaker
import pandas as pd
import sys
import os
from matplotlib import pyplot as plt
import numpy as np
import csv

plt.switch_backend('agg')
def main(argv = sys.argv):
    fig, (ax_q1, ax_i1)= plt.subplots(nrows=2,gridspec_kw={'height_ratios': [2, 2]})
    ofname = argv[3]+'/'+argv[1]+'_'+argv[2]+'.png'
    mfile1 = argv[1]+'_'+argv[2]+'_1bp_pwm.txt'
    mfile4 = argv[1]+'_'+argv[2]+'_mdsvm_pwm.txt'
    pfile = argv[1]+'_'+argv[2]+'_1bp_PWM_locs.out'
    fidp = open(pfile,'r')
    pwm=csv.reader(fidp,delimiter='\t')
    ploc = []
    pb = []
    pe = []
    pl = []
    ps = []
    for i,l,j,k,s in pwm:
        ploc.append(str(i))
        pl.append(int(l))
        pb.append(int(j))
        pe.append(int(k))
        ps.append(float(s))
    d = pd.read_csv(mfile1, delim_whitespace=True, index_col=0, comment='#')
    mat_q1 = logomaker.transform_matrix(pd.read_csv(mfile1, delim_whitespace=True, index_col=0, comment='#'),from_type="probability",to_type="information")
    logo_q1 = logomaker.Logo(mat_q1, ax = ax_q1, center_values=False, figsize = [100, 24])
    mat_i1 = logomaker.transform_matrix(pd.read_csv(mfile4, delim_whitespace=True, index_col=0, comment='#'))
    logo_i1 = logomaker.Logo(mat_i1, ax = ax_i1, center_values=False, figsize = [100, 24])
    for i in range(len(ploc)):
        if pb[i] > 0 and pe[i] <= 290:
            if i%2 == 1:
                ax_q1.text(pb[i], 2.3, ploc[i],fontsize=6)
            else:
                ax_q1.text(pb[i], 3, ploc[i],fontsize=6)
    ax_q1.set_ylim([0,4])
    ax_q1.set_ylabel("PWMs",fontsize=9)
    ax_i1.set_ylabel("mean dsvm",fontsize=9)
    fig.set_size_inches(40, 4)
    fig.tight_layout()
    fig.savefig(ofname)

main()
