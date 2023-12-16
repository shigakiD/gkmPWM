import logomaker
import pandas
import sys
import os
from matplotlib import pyplot as plt
import numpy as np
import csv
import argparse

plt.switch_backend('agg')
def main():
    parser = argparse.ArgumentParser(description='Plots of the TF binding profile of mapTF of a particular region')
    parser.add_argument('--fasta'  , required=True, type=str,
                        help='fasta file that was used when running mapTF')
    parser.add_argument('--locsprefix', required=True, type=str,
                        help='prefix of *_kmer_PWM_locs.out and *_motifs.out files')
    parser.add_argument('--weights', required=True, type=str,
                        help='kmer weight file that was used when running mapTF')
    parser.add_argument('--seqindex', required=True, type=str,
                        help='index of the sequence from the fasta file that will be plotted (first column in *_kmer_PWM_locs.out)')
    args = parser.parse_args()
    locsprefix = args.locsprefix
    bed = args.fasta
    sfile = args.fasta
    wfile = args.weights
    pfile = args.locsprefix+'_kmer_PWM_locs.out'
    mfile = args.locsprefix+'_motifs.out'
    ofname = args.locsprefix+'_'+args.seqindex+'_profile.png'
    seqnum = int(args.seqindex)
    fid = open(sfile,'r')
    for i in range(2*seqnum):
        seq = fid.readline().strip()
    fid.close()
    fig, (ax_s1,ax_q1, ax_i1)= plt.subplots(nrows=3,gridspec_kw={'height_ratios': [2, 2, 2]})
    fidp = open(pfile,'r')
    pwm=csv.reader(fidp,delimiter='\t')
    seqnumlist = []
    ploc = []
    pnum = []
    pb = []
    pe = []
    pd = []
    ps =[]
    sitelist = []
    for i,i2,i3,i4,i5,i6,i7,i8 in pwm:
        seqnumlist.append(int(i))
        ploc.append(str(i2))
        pnum.append(int(i3)-1)
        pb.append(int(i4))
        pe.append(int(i5))
        pd.append(float(i6))
        ps.append(float(i7))
    f = []
    i = 0
    lseq = len(seqnumlist)
    while seqnumlist[i] <= seqnum:
        if seqnumlist[i] == seqnum:
            f.append(i)
        i = i+1
        if i == lseq:
            break
    f.append(i)
    ploc = ploc[f[0]:f[-1]]
    pnum = pnum[f[0]:f[-1]]
    pb = pb[f[0]:f[-1]]
    pe = pe[f[0]:f[-1]]
    pd = pd[f[0]:f[-1]]
    ps = ps[f[0]:f[-1]]
    sitelist = sitelist[f[0]:f[-1]]
    names, motifs = getmotifs(mfile)
    fidp.close()
    [DSVM,l] = dsvm(wfile,seq) 
    [PWMmat,SEQ] = makePWMmatrix(pnum,pb,pe,motifs,seq)
    mat_s1 = logomaker.transform_matrix(pandas.DataFrame(SEQ, columns=['A', 'C', 'G', 'T']))
    logo_s1 = logomaker.Logo(mat_s1, ax = ax_s1, center_values=False, figsize = [100, 24])
    mat_q1 = logomaker.transform_matrix(pandas.DataFrame(PWMmat, columns=['A', 'C', 'G', 'T']),from_type="probability",to_type="information")
    
    logo_q1 = logomaker.Logo(mat_q1, ax = ax_q1, center_values=False, figsize = [100, 24])
    mat_i1 = logomaker.transform_matrix(pandas.DataFrame(DSVM, columns=['A', 'C', 'G', 'T']))
    logo_i1 = logomaker.Logo(mat_i1, ax = ax_i1, center_values=False, figsize = [100, 24])
    L = len(seq)
    for i in range(len(pb)):
        if pb[i] > 0 and pe[i] <= L-l+1:
            if i%2 == 1:
                ax_q1.text(pb[i], 2.3, ploc[i],fontsize=10)
            else:
                ax_q1.text(pb[i], 3, ploc[i],fontsize=10)
    ax_q1.set_ylim([0,4])
    ax_s1.set_ylabel("Sequence",fontsize=9)
    ax_q1.set_ylabel("PWMs",fontsize=9)
    ax_i1.set_ylabel("mean dsvm",fontsize=9)
    fig.set_size_inches(L/15, 6)
    fig.tight_layout()
    fig.savefig(ofname,format='png')

def getmotifs(mfile):
    fid = open(mfile, 'r')
    lines = fid.readlines()
    i = 0
    motifs = []
    names = []
    L = len(lines)
    while i < L:
        mat = []
        i = i+1
        names.append(lines[i].strip())
        i = i+1
        l = int(lines[i])
        i = i+1
        for i2 in range(l):
            m = lines[i].strip().split(' ')
            mrow = [float(j) for j in m]
            mat.append(mrow)
            i = i+1
        motifs.append(mat)
        i = i +1
    fid.close()
    return names, motifs

def makePWMmatrix(pnum,pb,pe,mat,seq):
    L = len(seq)
    PWM = []
    SEQ = []
    for i in range(L):
        PWM.append([0.25,0.25,0.25,0.25])
        vec = [0,0,0,0]
        if seq[i] == 'A':
            vec[0] = 1
        elif seq[i] == 'C':
            vec[1] = 1
        elif seq[i] == 'G':
            vec[2] = 1
        else:
            vec[3] = 1
        SEQ.append(vec)
    l = len(pnum)
    for i in range(l):
        l2 = pe[i]-pb[i]+1
        for j in range(l2):
            PWM[pb[i]+j-1][0] = mat[pnum[i]][j][0]
            PWM[pb[i]+j-1][1] = mat[pnum[i]][j][1]
            PWM[pb[i]+j-1][2] = mat[pnum[i]][j][2]
            PWM[pb[i]+j-1][3] = mat[pnum[i]][j][3]
            if seq[pb[i]+j-1] == 'A':
                SEQ[pb[i]+j-1][0] = 3
            elif seq[pb[i]+j-1] == 'C':
                SEQ[pb[i]+j-1][1] = 3
            elif seq[pb[i]+j-1] == 'G':
                SEQ[pb[i]+j-1][2] = 3
            else:
                SEQ[pb[i]+j-1][3] = 3

    return PWM,SEQ


def dsvm(wfile,seq):
    fid = open(wfile,'r')
    weights=csv.reader(fid,delimiter='\t')
    kmers = []
    scores = []
    for i,j in weights:
        kmers.append(str(i))
        scores.append(float(j))
    d = dict(zip(kmers,scores))
    l = len(kmers[1])
    L = len(seq)
    L2 = L-1
    S = []
    for i in range(L):
        S.append([0,0,0,0])
        if i-l < 0:
            loc = i
            start = 0
            end = l+i
            kmerlen = i
        elif i+l > L2:
            end = L
            start = i-l+1
            loc = l-1
            kmerlen = L-i-1
        else:
            start = i-l+1
            end = i+l
            loc = l-1
            kmerlen = l
        kmer = seq[start:end]
        ref = 0
        for j in range(kmerlen):
            kmer2 = kmer[(j):(l+j)]
            if kmer2 in d:
                ref = ref + d[kmer2]
            else:
                ref = ref + d[revcomp(kmer2)]
        kmerv = makevar(kmer,loc)
        alt = 0;
        for i2 in range(3):
            for j in range(kmerlen):
                kmer2 = kmerv[i2][(j):(l+j)]
                if kmer2 in d:
                    alt = alt + d[kmer2]
                else:
                    alt = alt + d[revcomp(kmer2)]
        if seq[i] == 'A':
            S[i][0] = ref-alt/3.0
        elif seq[i] == 'C':
            S[i][1] = ref-alt/3.0
        elif seq[i] == 'G':
            S[i][2] = ref-alt/3.0
        else:
            S[i][3]= ref-alt/3.0
    fid.close()
    return S,l
        
def revcomp(kmer):
    L = len(kmer)
    s=''
    for i in range(L):
        if kmer[i] == 'A':
            s='T'+s
        elif kmer[i] == 'C':
            s='G'+s
        elif kmer[i] == 'G':
            s='C'+s
        else:
            s='A'+s
    return s

def makevar(kmer,loc):
    l = len(kmer)
    kmerv=[]
    s = kmer
    loc2=loc+1
    if kmer[loc] == 'A':
        kmerv.append(s[:loc] + 'C' + s[loc2:])
        kmerv.append(s[:loc] + 'G' + s[loc2:])
        kmerv.append(s[:loc] + 'T' + s[loc2:])
    elif kmer[loc] == 'C':
        kmerv.append(s[:loc] + 'A' + s[loc2:])
        kmerv.append(s[:loc] + 'G' + s[loc2:])
        kmerv.append(s[:loc] + 'T' + s[loc2:])
    elif kmer[loc] == 'G':
        kmerv.append(s[:loc] + 'A' + s[loc2:])
        kmerv.append(s[:loc] + 'C' + s[loc2:])
        kmerv.append(s[:loc] + 'T' + s[loc2:])
    else:
        kmerv.append(s[:loc] + 'A' + s[loc2:])
        kmerv.append(s[:loc] + 'C' + s[loc2:])
        kmerv.append(s[:loc] + 'G' + s[loc2:])
        
    return kmerv

main()
