# gkmPWM

A package to extract compact and interpretable features from sequence based models.

## Installation

Run the following on the commandline:

```bash
git clone https://github.com/dshigaki/gkmPWM.git
```
This repository is organized into 4 directories.  Two of them, matlab and src contain the code to run:

<b>gkmPWMlasso</b>: an algorithm to extract known PWMs from a sequence based model.  

<b>gkmPWM</b>: an algorithm to learn <i>de novo</i> PWMs from a sequence based model.  

<b>mapTF</b>: a method to map the PWMs from gkmPWMlasso and gkmPWM to a set of sequences.  

To run the matlab code, include <i>addpath('dir/gkmPWM/matlab')</i> in one of your lines.  dir is the location of the gkmPWM directory.  These require MATLAB's statistics and machine learning toolbox.

To run the C code, in the gkmPWM directory, run
```bash
make 
```
There will be three executables: gkmPWMlasso, gkmPWM, and mapTF that appear in the directory.

There isn't really any advantage over running the MATLAB script over the C script and vice versa.  They run with comparable computation time and memory.  
  
In the visualization directory, there are two python scripts that create visuals for the outputs of the MATLAB and C functions.  Python3 is required to run these.  

Lastly, there is a directory containing example files that are used in the tutorial below.

The following sections will outline and provide examples of running gkmPWMlasso, gkmPWM, and mapTF.  Only required parameters and a handful of optional parameters will be discussed for the tutorial.  For information on the optional parameters, see the help sections in the code.  There will be lines of code that can be copied to generate the same outputs as in the example_files directory.  All file inputs in the examples are also in that directory.

## Running gkmPWMlasso

This function extracts motifs from a gkmSVM model from the R package or lsgkm (https://github.com/Dongwon-Lee/lsgkm).  It requires 3 parameters:
1. gkmSVM model prefix: <i>fileprefix</i>_svseq.fa and <i>fileprefix</i>_svalpha.out, or <i>fileprefix</i>.model.txt
2. Database of PWMs in meme format: <i>memefile</i>
3. number of motifs: _m_

The first input is the prefix of the output of gkmSVM or lsgkm.  The second input It uses a database of known motifs in meme format.  We provide a meme file, <i>combined_db_v4.meme</i>, containing nearly 2000 motifs.  The last input is the number of motifs to extract, which we will denote as _m_.  If you set this to 0, gkmPWMlasso will estimate the number of motifs to extract.  

<b>MATLAB</b>
```bash
gkmPWMlasso(fileprefix, memefile,m)
gkmPWMlasso('GM12878', 'combined_db_v4.meme', 30)
```
<b>C</b>
```bash
./gkmPWMlasso fileprefix memefile m
./gkmPWMlasso GM12878 combined_db_v4.meme 30
```

These both output <i>GM12878_10_6_30_gkmPWMlasso.out</i>, which contains the following columns:
1. <u>Cluster ID</u>: The cluster to which the PWM belongs.  gkmPWMlasso clusters PWMs to prevent linear dependence and redundant features.  
2. <u>ID</u>: The number of the motif from (1) as it as appears in the memefile input.
3. <u>MOTIF ID</u>: The name of the motif in the memefile
4. <u>Weight</u>: The regression weight of the motif.
5. <u>Z-score</u>: The number of standard deviations of the average of the top gapped k-mers weights.  The number of gapped k-mers used depends on the combination of optional parameters '<i>l</i>', '<i>k</i>', and '<i>KmerFrac</i>'.
6. <u>Importance</u>: The relative increase in error when removing that PWM from the list of features.

You can create a pdf of the output by running <i>plotMotif.py</i> in the visualization directory.  The parameters required are:
1. gkmPWMlasso output: <i>lassofile</i> passed to --info
2. Database of PWMs in meme format: <i>memefile</i> passed to --meme
3. Output prefix: _outprefix_ passed to --output
```bash
python plotMotif.py --info lassofile --meme memefile --output outprefix
python plotMotif.py --info GM12878_10_6_30_gkmPWMlasso.out --meme combined_db_v4.meme --output GM12878_10_6_30_gkmPWMlasso
```
This will create _GM12878_10_6_30_gkmPWMlasso.pdf_
![](https://github.com/shigakiD/gkmPWM/blob/main/example_files/GM12878_10_6_30_gkmPWMlasso.png)
A quick note:  for the optional parameters <i>l</i> and <i>k</i>, these <b>do not</b> need to be the same as the <i>l</i> and <i>k</i> from the gkmSVM model.  In fact, if <i>l</i> and <i>k</i> generate too many gapped k-mers, gkmPWMlasso will take only a subset of the gapped k-mer to use as features.  This is also true for gkmPWM.

## Running gkmPWM (de novo PWMs)

This function learns PWMs <i>de novo</i> from sequence based models.  The required parameters are:  
1. gkmSVM model prefix: <i>fileprefix</i>_svseq.fa and <i>fileprefix</i>_svalpha.out, or <i>fileprefix</i>.model.txt
2. gkmSVM kmer weights: <i>wfilename</i> (this is used to seed the motifs to make it converge faster) 
3. Database of PWMs in meme format: <i>memefile</i>
4. number of motifs: <i>m</i> 

Inputs 1, 3, and 4 are the same as the inputs to gkmPWMlasso.  

Input 2 is a kmer weight file that can be generated using the code from lsgkm.  You can generate these using the following commands.  
```bash
python /lsgkm/scripts/nrkmers.py l outputname
```
'<i>l</i>' is the length of the kmer of interest.  outputname is a fasta file containing all length <i>l</i> kmers.  Reverse complements are treated as similar.  To generate the weight file, run:
```bash
/lsgkm/src/gkmpredict outputname fileprefix.model.txt wfilename
```
The output is a list of scores for each full length kmer.

A weight file (GM12878_weights.out) is provided for you.  You can generate this yourself by using <i>l=11</i> and using GM12878.model.txt as the svm model  

To run gkmPWM:

<b>MATLAB</b>
```bash
gkmPWM(fileprefix, wfilename, memefile, m)
gkmPWM('GM12878', 'GM12878_weights.out', 'combined_db_v4.meme',15)
```
<b>C</b>
```bash
./gkmPWM fileprefix wfilename memefile m
./gkmPWM GM12878 GM12878_weights.out combined_db_v4.meme 15
```
This creates 3 files <i>GM12878_10_6_0_15_denovo.meme</i>,<i>GM12878_10_6_0_error.out</i>, and  <i>GM12878_10_6_0_15_gkmPWM.out</i>.  

The first is a meme file containing the PWMs.  The second file is a summary file.  The second file is a record of the error after each iteration.  If the error does not clearly converge, run it again for more reps.  There will be small jumps in the error, which is gkmPWM adjusting the alignment of the PWMs.   The last file is a summary file with the following columns:

1. <u>MOTIF</u>: The name of the motif in the memefile that had the highest pearson correlation with that PWM
2. <u>ID</u>: The number of the motif from (1) as it as appears in the memefile input.
3. <u>Similarity</u>: The pearson correlation from (1).
4. <u>Redundancy</u>: The highest pearson correlation with another PWM that was also learned.
5. <u>Weight</u>: The regression weight of the motif.
6. <u>Z-score</u>: The number of standard deviations of the average of the top gapped k-mers weights.  The number of gapped k-mers used depends on the combination of optional parameters '<i>l</i>', '<i>k</i>', and '<i>KmerFrac</i>'.
7. <u>Importance</u>: The relative increase in error when removing that PWM from the list of features.

The numbers in the output correspond to the '<i>l</i>', '<i>k</i>', and '<i>RegFrac</i>' parameters.

If you want to increase the information in the PWMs (make them less noisy), you can change the optional parameter '<i>RegFrac</i>' in MATLAB and '<i>-r</i>' in C.  It can take values in [0,1).  The default value is 0.  From my own experimentation, this value works best in the interval [0, 0.05].

You can create a pdf of the output by running <i>plotMotif.py</i> in the visualization directory.  The parameters required are:
1. gkmPWM information output: <i>gkmPWMfile</i> passed to --info
2. gkmPWM meme output: <i>memefile</i> passed to --meme
3. Output prefix: _outprefix_ passed to --output
You also need to specify --denovo
```bash
python plotMotif.py --denovo --info gkmPWMfile --meme memefile --output outprefix
python plotMotif.py --denovo --info GM12878_10_6_0_15_gkmPWM.out --meme GM12878_10_6_0_15_denovo.meme --output GM12878_10_6_0_15_gkmPWM
 ```
 This will create _GM12878_10_6_0_15_gkmPWM.pdf_
## Running mapTF 

This function maps PWMS from both gkmPWM and gkmPWMlasso to sequences at base-pair resolution.  The required parameters are:
1. Sequences in fasta format: <i>seqfile</i>.
2. gkmSVM kmer weights: <i>wfilename</i> 
3. denovo.meme output of gkmPWM: <i>denovofile</i>
4. gkmPWMlasso out: <i>lassofile</i>
5.  Database of PWMs in meme format: <i>memefile</i>
6. Output prefix: <i>outprefix</i>

 <b>The optional parameters '<i>l</i>' and '<i>k</i>' must be the same as the model's parameters</b>.
 
 To run mapTF:
 
<b>MATLAB</b>
```bash
mapTF(seqfile, wfilename, denovofile, lassofile, memefile,outprefix)
gkmPWM('GM12878.fa', 'GM12878_weights.out','GM12878_10_6_0_15_denovo.meme', 'GM12878_10_6_30_gkmPWMlasso.out', 'combined_db_v4.meme','GM12878')
```
<b>C</b>
```bash
./mapTF seqfile wfilename denovofile lassofile memefile outprefix
./mapTF GM12878.fa GM12878_weights.out GM12878_10_6_0_15_denovo.meme GM12878_10_6_30_gkmPWMlasso.out combined_db_v4.meme GM12878
```
This outputs two files:
<i>outprefix_motifs.out</i>
<i>outprefix_kmer_PWM_locs.out</i>
The first file contains the PWMs for the visualization script mapTF_profile.py.  
The second file gives a list of the locations of the mapped TFBSs.  The columns are:

1. <u>Sequence ID</u>: The sequence index to which the TFBS was mapped (not zero indexed)
2. <u>Motif Name</u>
3. <u>Motif number</u> in <i>outprefix_motifs.out</i>   
4. <u>Start location</u> (not zero indexed)
5. <u>End location</u>
6. <u>Average kmer probability</u> (higher means more likely a binding site)
7. <u>Correlation with deltaSVM</u> (higher means more likely a good match).
8. <u>TFBS sequence</u>

You can convert the _outprefix_kmer_PWM_locs.out_ file to a bed by using the _convert2bed.m_ or <i>convert2bed.py</i> functions.  The _bedfile_ input should be in the same order as the input fasta file.

<b>MATLAB</b>
```bash
convert2bed(outprefix,bedfile)
convert2bed('GM12878', 'GM12878.bed')
```
<b>Python</b>
```bash
python convert2bed.py --locsprefix outprefix --bed bedfile
python convert2bed.py --locsprefix GM12878 --bed GM12878.bed
```
You can make a profile plot of the results of mapTF for a sequence using <i>mapTF_profile.py</i> for a given sequence.  It requires: 

1. Sequences in fasta format: <i>seqfile</i>.
2. gkmSVM kmer weights: <i>wfilename</i> 
3. mapTF output prefix: <i>kmerPWMprefix</i>
4. Sequence Index: _sind_ (from _kmer_PWM_locs.out_)
```bash
python mapTF_profile.py --fasta seqfile --weights wfilename --locsprefix kmerPWMprefix --seqindex sind
python mapTF_profile.py --fasta GM12878.fa --weights GM12878_weights.out --locsprefix GM12878 --seqindex 1
```
This creates a png named _GM12878_1_profile.png_, formatted as _kmerPWMprefix_sind_profile.png_.  The first row is the sequences of interest with the nucleotides contained in binding sites raised.  The second is are the mapped PWMs.  The third row is the average deltaSVM score of each nucleotide for all possible mutations.

## Authors

* **Dustin Shigaki** 
* **Gary Yang** 
* **Michael A Beer** 
