import sys, os, getopt

def main():
    f = ''
    w = ''
    n =''
    m =''
    l='10'
    k='6'
    r='0'
    c='true'
    argv = sys.argv[1:]
    opts, arg = getopt.getopt(argv,'f:w:n:m:k:l:r:c:')
    for opt, arg in opts:
        if opt in ("-f"):
            f = arg
        elif opt in ("-w"):
            w = arg
        elif opt in ("-n"):
            n = arg
        elif opt in ("-m"):
            m = arg
        elif opt in ("-l"):
            l = arg
        elif opt in ("-k"):
            k = arg
        elif opt in ("-r"):
            r = arg
        elif opt in ("-c"):
            c = arg
    ofile = open('qgkmPWM','w')
    ofile.write('#!/bin/bash\n')
    ofile.write('#SBATCH --time=72:0:0\n')
    ofile.write('#SBATCH --mem=10G\n')
    ofile.write("matlab -nojvm -singleCompThread -nodesktop -nosplash -r "+'"'+"gkmPWM('"+f+"','"+w+"',"+m+","+n+",'l',"+l+",'k',"+k+","+"'RegFrac',"+r+",'RC',"+c+")"+'"'+"\n")
    ofile.close()
    os.system("sbatch qgkmPWM")
main()
