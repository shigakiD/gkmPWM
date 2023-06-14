import sys, os

def main(argv=sys.argv):
    ofile = open('qLASSO','w')
    ofile.write('#!/bin/bash\n')
    ofile.write('#SBATCH --time=2:0:0\n')
    ofile.write('#SBATCH --mem=12G\n')
    ofile.write("matlab -nojvm -nodesktop -nosplash -r "+'"'+"tic;gkmPWMlasso2('"+argv[1]+"',"+argv[2]+",'"+argv[3]+"');toc"+'"'+"\n")
    #ofile.write("matlab -nodesktop -nosplash -r "+'"'+"plot_motif_lasso2('"+argv[1]+"',"+argv[2]+");quit"+'"'+"\n")
    ofile.close()
    os.system("sbatch qLASSO")
main()
