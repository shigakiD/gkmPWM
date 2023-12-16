import sys,os,argparse

def main(argv = sys.argv):
    parser = argparse.ArgumentParser(description='Converts the results of mapTF *_kmer_PWM_locs.out to bed format')
    parser.add_argument('--bed'  , required=True, type=str,
                        help='bed file in the same order of the fasta file input to mapTF')
    parser.add_argument('--locsprefix', required=True, type=str,
                        help='prefix of *_kmer_PWM_locs.out file.')
    args = parser.parse_args()
    locsprefix = args.locsprefix
    bed = args.bed
    fid = open(locsprefix+'_kmer_PWM_locs.out', 'r')
    x = fid.readlines()
    fid.close()
    fid = open(bed, 'r')
    y = fid.readlines()
    fid.close()
    L = len(x)
    fid = open(locsprefix+'_kmer_PWM_locs.bed', 'w')
    L2 = len(y)
    Y = []
    for i in range(L2):
        Y.append(y[i].strip().split('\t'))
        Y[i][1] = int(Y[i][1])
        Y[i][2] = int(Y[i][2])
    for i in range(L):
        line = x[i].split('\t')
        n = int(line[0])-1
        b = int(line[3])+Y[n][1]-1
        e = int(line[4])+Y[n][1]
        fid.write(Y[n][0]+'\t'+str(b)+'\t'+str(e)+'\t'+line[1]+'\t'+line[5]+'\t'+line[6]+'\t'+line[7])
    fid.close()

main()
