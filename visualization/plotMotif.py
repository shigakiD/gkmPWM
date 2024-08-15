import sys, re, io, argparse
from typing import *
import pandas as pd
import numpy as np
import logomaker as lm

from PIL import Image

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib
matplotlib.use('Agg')

def parsePWMsFromMeme(fname: str) -> Dict[int, Tuple[str, np.ndarray]]:
    """ This function takes in a meme format file and
    outputs a dictionary with key being 1-based counter
    and value being a tuple (motif name, motif PWM [np.ndarray])
    """
    
    def getMotif(lines: List[str]) -> Tuple[str, np.ndarray]:
        """ Files in the .meme format contain many motif PWMs.
            Besides headers, each moti is described by the following
            lines:
                1. Motif name
                2. Description, including PWM length 'w' and PWM alphabet size 'l'
                3. PWM itselfs (spanning 'w' lines)
                4. Empty line
        """
        motifName = lines[0].strip().split(' ')[-1]
        pwm = lines[2:]
        pwm = [[float(m) for m in re.split('  | |\t', position.strip())] 
               for position in pwm 
               if position.strip()]
        motifLength = re.findall('w= [0-9]+', lines[1])[0][3:]
        assert len(pwm) == int(motifLength)
        return (motifName, np.array(pwm))    
    
    with open(fname) as fp:
        mfile = fp.readlines()
    motifNameLines = [idx for idx, line in enumerate(mfile) 
                      if 'MOTIF' in line]
    motifPWMs = dict()
    for counter in range(len(motifNameLines[:-1])):
        begin = motifNameLines[counter]
        end = motifNameLines[counter+1]
        currMotif = mfile[begin:end]
        motifPWMs[counter+1] = getMotif(currMotif)
    motifPWMs[counter+2] = getMotif(mfile[motifNameLines[counter+1]:])
    return motifPWMs


def parseLassoOut(fname: str) -> pd.DataFrame:
    """ This function parses *_gkmPWMlasso.out 
    and ensures output dataframe contains 1 representative
    motif entry for each cluster of similar motifs """
    df = pd.read_csv(fname, 
                     sep="\t", 
                     skiprows=4, 
                     header=None, 
                     names=["clus", "uid", "motif", "W", "Z", "I"])
    uniqMotif = df.loc[np.concatenate([[1], np.diff(df.clus)]).astype(bool)]
    uniqMotif = uniqMotif.sort_values(by="Z", ascending=False)
    return uniqMotif


def parseDeNovoOut(fname: str) -> pd.DataFrame:
    """ This function parses *_gkmPWM.out 
    and ensures output dataframe contains a 1-based 
    unique ID column 
    """
    df = pd.read_csv(fname, 
                     sep="\t", 
                     skiprows=3, 
                     header=None, 
                     names=["motif", "uid", "R", "NA", "W", "Z", "I"])
    df.uid = np.arange(1, df.shape[0]+1)
    return df


def generate_pdf(out: pd.DataFrame, 
                 pwms: Dict[int, Tuple[str, np.ndarray]], 
                 prefix: str, 
                 denovo=False):
    
    numMotif = out.shape[0]
    ptScale  = 85
    offsets  = 0.25
    headerY  = (numMotif + offsets) * ptScale
    
    """ Set up canvas """
    fig = plt.figure(figsize=(8, 7/30*numMotif), dpi=200)
    cax = plt.gca()
    plt.xlim(-1200, 1600)
    plt.ylim(- 100, (numMotif+1)*ptScale)
    
    """ Visualize header """
    def writeHeader(x, text, y=headerY, ax=cax, s=ptScale):
        ax.text(x*s, y, text)
    
    """ For each motif entry, visualize
            1. Heatmap for quantitative info (regression weight, z-score, error),
            2. Quantitative information
            3. Non-quantitative information (motif name, unique ID)
            4. PWM logo
    """
    def plotHeatmap(col, idx, x, y, df=out, ax=cax, s=ptScale):
        value = df[col].iloc[idx]
        maxV  = df[col].abs().max()
        color = ([1-value/maxV, 1-value/maxV, 1] 
                 if value > 0 
                 else [1, 1+value/maxV, 1+value/maxV])
        ax.add_patch(
            Rectangle((x*s, y*s), s, s, fc=color, ec='k', ls='-', lw=0.5)
        )
        
    def writeQInfo(x, y, col, df=out, ax=cax, s=ptScale):
        ax.text(x*s, y*s, f"{df[col].iloc[idx]:>5.2f}", size=9)
        
    def writeInfo (x, y, text, df=out, ax=cax, s=ptScale):
        ax.text(x*s, y*s, text, size=9)

    def vizPWM(pwm, x, y, s=ptScale, factor=100, ax=cax):
        pwm_info = lm.transform_matrix(
            pd.DataFrame(pwm/pwm.sum(axis=1, keepdims=True), 
                         columns=['A', 'C', 'G', 'T']),
            from_type="probability",
            to_type="information")
        figTmp, axTmp = plt.subplots(1, 1, figsize=(0.4*pwm.shape[0], 1), dpi=200)
        lmObj = lm.Logo(pwm_info, ax=axTmp, center_values=False)
        lmObj.draw()
        axTmp.tick_params(bottom=False, labelbottom=False) 
        axTmp.spines[  'top'].set_visible(False)
        axTmp.spines['right'].set_visible(False)
        figTmp.tight_layout()
        b = io.BytesIO()
        figTmp.savefig(b, format='png', dpi=200, bbox_inches="tight", pad_inches=0)
        b.seek(0)
        img = Image.open(b)
        ax.imshow(img, extent=[x*s, (x+0.4*pwm.shape[0])*s, y*s, (y+1)*s], 
                  clip_on=True, interpolation="antialiased", interpolation_stage="data")  
        plt.close(figTmp)
    

    writeHeader(-11.6 if denovo else -9   , 'Motif')
    writeHeader(- 4.3 if denovo else -11.6, 'R' if denovo else 'ID')
    for x, t in zip([-12.8, -1.8, 0.6, 3.2, 5.5, 6.6, 7.7], 
                    ['N', 'W', 'Z', 'I', 'W', 'Z', 'I']):
        writeHeader(x,t)
        
    for idx in range(numMotif):
        for col, x in zip(['W', 'Z', 'I'], [5.3, 6.3, 7.3]):   
            plotHeatmap(col, idx, x, numMotif-idx-1)
            
        motifUID = out.iloc[idx].uid
        infoYPos = numMotif-idx-0.8
        
        (name, pwm) = pwms[motifUID]
        name = name[:min(13, len(name))]
        writeInfo(-11.6 if denovo else -9, infoYPos, name)
        writeInfo(-12.8, infoYPos, idx+1)

        for col, x in zip(['W', 'Z', 'I'], [-2, 0.4, 2.6]):
            writeQInfo(x, infoYPos, col)
        if denovo:
            writeQInfo(-4.5, infoYPos, 'R')
        else:
            writeInfo(-11.6, infoYPos, motifUID)
        
        vizPWM(pwm, 9.1, numMotif-idx-1.25)
        
    plt.axis('off')
    try:
        plt.savefig(f"{prefix}.pdf" if len(prefix) > 4 and prefix[-4:] != '.pdf' else prefix, dpi=600)
    except:
        print("ERROR: Fail to write a pdf to the current directory.")


def main():
    parser = argparse.ArgumentParser(description='Visualize PWM logos and regression metrics for gkmPWMlasso and gkmPWM')
    parser.add_argument('--denovo', default=False, action='store_true', 
                        help="specify this, if visualizing gkmPWM's output")
    parser.add_argument('--meme'  , required=True, type=str,
                        help='filename to a meme file. If visualizing gkmPWM output, this should be *_denovo.meme' + 
                        'If visualizing gkmPWMlasso output, this should be the motif database file used during training.')    
    parser.add_argument('--info'  , required=True, type=str,
                        help='filename to the output of either gkmPWMlasso (*_gkmPWMlasso.out) or gkmPWM (*_gkmPWM.out)')
    parser.add_argument('--output', required=True, type=str,
                        help='prefix to generated pdf file.')
    args = parser.parse_args()
    
    if args.denovo:
        if '_gkmPWM.out' not in args.info or '_denovo.meme' not in args.meme:
            print(f"ERROR: Since you selected to visualize output of gkmPWM by specifying \n" + 
                  f"       --denovo, your meme file should end in _denovo.meme and your output\n" + 
                  f"       file should end in _gkmPWM.out. Currently they are \n" + 
                  f"       {args.meme} \n" + 
                  f"       and \n" + 
                  f"       {args.info}")
            return
    else:
        if '_gkmPWMlasso.out' not in args.info:
            print(f"ERROR: Since you selected to visualize output of gkmPWMlasso by NOT specifying \n" + 
                  f"       --denovo, your output file should end in _gkmPWMlasso.out. Currently, it \n" + 
                  f"       is {args.info}. \n" + 
                  f"       Also, please specify the motif database file you used for gkmPWMlasso \n" + 
                  f"       with the option -meme.")
            return
    
    print("Visualizing gkmPWM de novo motifs" if args.denovo else "Visualizing gkmPWMlasso motifs")
    
    try:
        memeFile = parsePWMsFromMeme(args.meme)
    except:
        print(f"ERROR: The meme file is not properly formatted")
        print(f"The meme file should provide motif PWMs. Besides headers, each motif is described by the following lines: \n" + 
              f"      1. Motif name\n" +
              f"      2. Description, including PWM length 'w' and PWM alphabet size 'l'\n" +
              f"      3. PWM itselfs (spanning 'w' lines)\n" +
              f"      4. An empty line\n" + 
              f"Please see the provided combined_db_v4.meme for examples.")
        return
        
    try:
        outpFile = parseDeNovoOut(args.info) if args.denovo else parseLassoOut(args.info)
    except:
        print(f"ERROR: The output file is corrupted")
        return
    
    generate_pdf(outpFile, memeFile, args.output, args.denovo)
    

if __name__ == '__main__':
    main()
