# gkmPWM

A method to extract compact and interpretable features from a gkmSVM model.  The outputs are PWMs in the meme format.

## Getting Started

Run the following on the commandline:

```bash
git clone https://github.com/dshigaki/gkmPWM.git
```

Requires MATLAB's statistics and machine learning toolbox.  Some python code are included.  You can go through the examples below and compare to the example outputs given (the files with "example").

I highly recommend using the default settings for each function, especially for the parameters "l" and "k".  There is a trade-off between the number of gapped k-mers features and the amount of memory and computation time needed.  I find that (l,k) = (10,6) provides the ideal balance.

## Running gkmPWMlasso

This function requires either a gkmSVM model from the R version, <i>filename_svseq.fa</i> and <i>filename_svalpha.out</i>, or the lsgkm version, <i>filename.model.txt</i>.  The following line finds _m_ motifs a given _memefile_:

```bash
python run_gkmPWMlasso.py filename m memefile 
```

A memefile <i>combined_db_v4.meme</i> is provided containing approximately 2000 motifs. Some of them are quite similar, but gkmPWMlasso clusters the similar PWMs together to avoid linear dependence.

```bash
python run_gkmPWMlasso.py GM12878 30 combined_db_v4.meme
```
This gives us <i>GM12878_10_6_30_gkmPWMlasso.out</i>  
  
If you want to create a pdf with the motif logos on them, then run 

```bash
python run_gkmPWMlasso_pdf.py filename m memefile  
python run_gkmPWMlasso_pdf.py GM12878 30 combined_db_v4.meme
```

This will create a pdf of the <i>motifs.out</i> file, with only the logo of the motif at the top of the cluster.

This version is the faster of the two PWM methods, since it uses a set of PWMs from other databases.  It takes about 15-20 minutes for (l,k) = (10,6).  Also, I usually set the number of motifs to 20 or 30.  based

## Running gkmPWM (de novo motif learning)

This particular function requires:  
1. gkmSVM models: <i>filename_svseq.fa</i> and <i>filename_svalpha.out</i>, or <i>filename.model.txt</i>
2. gkmSVM kmer weights: <i>wfilename</i> (this is used to seed the motifs to make it converge faster) 
3. number of motifs: <i>m</i> 
4. number of iterations to run: <i>n</i> ( I recommend at least ten times the number of motifs)  

A weight file (GM12878_weights.out) is provided for you.  You learn out to generate these by going to the Beer lab website or Dongwon's github.  
Enter the following into the commandline:

```bash
python run_gkmPWM.py -f filename -w wfilename -m motif# -n iteration#
python run_gkmPWM.py -f GM12878 -w GM12878_weights.out -m 15 -n 200
```

This creates four files <i>GM12878_10_6_0_15_denovo.meme</i>, <i>GM12878_all_10_6_0_15_denovo.meme</i>, <i>GM12878_10_6_0_15_gkmPWM.out</i>, and <i>GM12878_10_6_0_error.out</i>.  The first is a meme file with a list of the most predictive motifs.  The second file is the full list of motifs.  The third file is a summary file.  The last file is a record of the error after each iteration.  If the error does not clearly converge, run it again for more reps.  There should also be little jumps in the error after 20 iterations.  This is the motif lengths being adjusted periodically, so don't worry.  

If you want to reduce the noise in the PWMs, you can run 

```bash
python run_gkmPWM.py -f filename -w wfilename -m motif# -n iteration# -r Rf
python run_gkmPWM.py -f GM12878 -w GM12878_weights.out -m 15 -n 200 -r 0.03
```

Rf is the degree of noise reduction.  It can take values in [0,1).  From my own experimentation, a good Rf to use lies between [0.005, 0.03].  

This creates four files <i>GM12878_10_6_0.03_15_denovo.meme</i>, <i>GM12878_all_10_6_0.03_15_denovo.meme</i>, <i>GM12878_10_6_0.03_15_gkmPWM.out</i>, and <i>GM12878_10_6_0.03_error.out</i>. 

## Mapping motifs to sequence

Before starting:
1.	Throw all the sequences that you want mapped into one fasta file.  
2.	Run gkmPWM denovo and gkmPWMlasso on the revelant dataset

NOTE: “ofn” should be the same for all functions  
  
**function seq2prob(sfn, wfn, ofn)**  
sfn: fasta file  
wfn: kmer weight file  
ofn: output header  
  
There is a directory fasta file called GM12878.fa.  It has only 5 sequences.  
Running the following line in MATLAB creates 5 files: GM12878_1_prob.out, GM12878_2_prob.out…

```bash
seq2prob('GM12878.fa', 'GM12878_weights.out', 'GM12878')
```

I also create variant scores.  
  
**function seq2var(sfn, wfn, ofn)**  
sfn: fasta file  
wfn: kmer weight file  
ofn: output header  
  
It has the same input format as the previous function, so just run the line above with the different function name  
  
This will generate GM12878_1_dsvm.out, GM12878_2_dsvm.out…  
  
These scores then need to converted into the proper format.  
  
**function avg_dsvm(fn,ofn)**  
fn:fasta file  
ofn: output fileheader  
  
Reminder: ofn should be the same as your input for seq2var.m.  

```bash
avg_dsvm('GM12878.fa', 'GM12878')
```

The result is GM12878_1_mdsvm_pwm.txt etc.  
  
Then, we need to select motifs.  This requires a motif file from gkmPWM and gkmPWMlasso.  
  
<b>function process_motifs(dfn, lfn, ofn)</b>  
dfn: file name for denovo motifs  
lfn: file name for lasso motifs  
ofn: output filename  
  
There is a file called GM12878_all.out which is the output of this function.  

We can now map the motifs to each sequence  

**function mapTF(fn, mfn, rnum, GC, l,k,ofn)**  
fn: fasta file used for seq2prob  
mfn: file name for motifs from process motifs  
rnum: sequence number in fasta file  
GC: GC content as a fraction between 0 to 1
l: kmer length
k: number of ungapped positions
ofn: output file header  
  
Entering the following line yields a few files which are needed to generate the nice plots  

```bash
mapTF('GM12878.fa', 'GM12878_all.out', 1, 0.46, 11,7,'GM12878')
```

GM12878_1_kmer_pwm.txt  
GM12878_1_kmer_PWM_locs.out  
  
The second file contains the locations of the binding sites  

I also created a version that is basically a PWM scan that does not require you to use kmer weights.  

**function mapTF_noweights(fn, mfn, r, rnum, GC, l,k,ofn)**
fn: fasta file used for seq2prob
mfn: file name for motifs from process motifs
r: the fraction of bases that are in a TFBS from 0 to 1
rnum: sequence number in fasta file
GC: GC content as a fraction between 0 to 1
l: kmer length
k: number of ungapped positions
There is an additional parameter "r" needed to run this function.  It's basically a number that controls how sensitive you want the motif mapping to be.  The larger r is, the more motifs that you will call.   

```bash
mapTF_noweights('GM12878.fa', 'GM12878_all.out', 0.1, 1, 0.46, 11,7,'GM12878')
```

GM12878_1_1bp_pwm.txt
GM12878_1_1bp_PWM_locs.out

The second file contains the locations of the binding sites


To generate the images,  

```bash
python mat2logo_kmer.py ofn rnum outdirectory
```

This script requires the logomaker and pandas packages in python.  

I usually create a new directory for all of the images to convenience.    

```bash
python mat2logo_kmer.py GM12878 1 images/ 
```

This creates images/GM12878_1.png  

## Authors

* **Dustin Shigaki** * 
* **Michael A Beer** *
