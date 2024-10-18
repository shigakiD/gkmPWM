//
//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
//
//
// Include files
#include "main.h"
#include "gkmPWM.h"
#include "gkmPWM_terminate.h"
#include "gkmPWM_emxAPI.h"
#include "gkmPWM_emxutil.h"
#include "cblas.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>  // for strtol
#include <ctype.h>
#include <errno.h>
#include <unistd.h>

emxArray_char_T* allocate_for_charArray(char* str) {
  emxArray_char_T *s1;
  emxInit_char_T(&s1, 2);
  int length = strlen(str);
  s1->data = (char*) malloc(sizeof(char)*length);
  memcpy(s1->data, str, length);
  s1->size[0] = 1;
  s1->size[1] = length;
  s1->allocatedSize = length;
  s1->canFreeData = true;
  return s1;
}
  
void display_arguments() {
    printf(
            "\n"
            "gkmPWM: A method to extract compact and predictive motifs from sequence-based models of regulatory elements de novo.\n"
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki, Gary Yang, Michael Beer\n"
            "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     gkmPWM [options] <prefix> <weight file> <database> <number of positive motifs>\n"
            "\n"
            "Arguments:\n"
            "                     prefix: prefix of a gkmSVM model, where the model files are either\n"
            "                             FILEHEADER_svseq.fa and FILEHEADER_svalpha.out OR FILEHEADER.model.txt\n"
            "                weight file: file name to the two-column, tab-delimited kmer weight file\n"
            "                   database: The collection of position weight matrix (PWM) in meme format. \n"
            "                             An example file (combined_db_v4.meme) is provided.\n"
            "  number of positive motifs: The number of positive set motifs to learn\n"
            "\n"
            "  Options:\n"
            "    -c <float>    If a PWM has a pearson correlation greater than this value with another PWM, it \n"
            "                  will be reseeded. This is to prevent linear dependence and redundancy. (default: 0.9)\n"
            "    -r <float>    This parameter will push the PWMs to have higher information. Must be a value in [0,1). \n"
            "                  Keeping this value low (<0.05) works best. (default: 0).\n"
            "    -P <float>    The ratio of positive PWMs to negative PWMs that will be seeded.  e.g., PNratio = 2 \n"
            "                  means twice as many positive set PWMs will be seeded as negative set PWMs. \n"
            "                  When set to 0, gkmPWM will automatically infer this number. (default: 0)\n"
            "    -f <float>    Set the fraction of the total number of gapped k-mers to use with gkmPWM. \n"
            "                  This reduces the memory and runtime needed. If the total number of gapped k-mers \n"
            "                  is too high with the given combination of (l,k,KmerFrac), KmerFrac will be \n"
            "                  automatically set to a lower value to create a more workable number of gapped k-mers \n"
            "    -l <int>      The full length of the gapped k-mer.  This DOES NOT need to be the same as the \n"
            "                  l in the gkmSVM model (default: 11) \n"
            "    -k <int>      The number of ungapped positions of the gapped k-mer. This DOES NOT need to be \n"
            "                  the same as the k in the gkmSVM model (default: 7)\n"
            "    -i <int>      Maximum number of iterations gkmPWM will run before exiting. If the loss \n"
            "                  converges, gkmPWM will finish prior to this (default: 200)\n"
            "    -B            if set, use both the positive and negative set to get the background \n"
            "                  distribution of gapped k-mers. This is best used when the model is trained \n"
            "                  with both the positive and negative sets containing regulatory elements.\n"
            "    -R            if set, consider reverse-complements of gapped k-mers to be distinct \n"
            "    -L            if set, WILL NOT automatically lower the fraction of gapped-kmer used. When set, \n"
            "                  the program will require more memory and need longer runtime \n"
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
    if(argc == 1) { display_arguments(); }

    double rCorr = 0.90;
    double reg = 0;
    double iPNRatio = 2;
    double lSVM = 11;
    double kSVM = 7 ;
    double numIterations = 200;
    double kmerFrac = 1;
    int backgroundGC = 0;
    int reverseCompl = 1;
    int kmerFracLimit = 1;
    char * pEnd;

    int c;
    while ((c = getopt(argc, argv, "BRLc:r:P:l:k:i:f:")) != -1) {
        switch (c) {
            case 'B':
                backgroundGC = 1;
                break;
            case 'R':
                reverseCompl = 0;
                break;
            case 'L':
                kmerFracLimit = 0;
                break;
            case 'c':
                rCorr = strtod(optarg, &pEnd);
                if (rCorr <= 0) {
                    printf("ERROR: Correlation cutoff must be a positive fraction.\n");
                    exit(1);
                }
                break;            
            case 'r':
                reg = strtod(optarg, &pEnd);
                if (reg >= 1 || reg < 0) {
                    printf("ERROR: regularization must be a fraction greater than or equal to 0 and less than 1\n");
                    exit(1);
                } else if (reg > 0.05) {
                    printf("Warning: Recommended values for RegFrac are in the interval [0, 0.05]. Still running with RegFrac = %f.\n", reg);
                }
                
                break;
            case 'P':
                iPNRatio = strtod(optarg, &pEnd);
                if (iPNRatio < 0) {
                    printf("ERROR: Positive-Negative PWM ratio must be a positive fraction.\n");
                    exit(1);
                }
                break;
            case 'l':
                lSVM = strtod(optarg, &pEnd);
                if (lSVM <= 0) {
                    printf("ERROR: l must be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'k':
                kSVM = strtod(optarg, &pEnd);
                if (kSVM <= 0) {
                    printf("ERROR: k must be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'i':
                numIterations = strtod(optarg, &pEnd);
                if (numIterations <= 0) {
                    printf("ERROR: Number of optimization iteration must be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'f':
                kmerFrac = strtod(optarg, &pEnd);
                if (kmerFrac <= 0 || kmerFrac > 1) {
                    printf("ERROR: Fraction K-mer used be a positive fraction in (0, 1].\n");
                    exit(1);
                }
                break;
            default:
                fprintf(stderr, "Unknown option: -%c\n", c);
                display_arguments();
        }
    }

    if (argc - optind != 4) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the documentation.\n");
        display_arguments();
    }
    if (lSVM < kSVM) {
        printf("ERROR: k must be a positive integer smaller than l\n");
        exit(1);
    }
    
    int index = optind;
    emxArray_char_T *model_file  = allocate_for_charArray(argv[index++]);
    emxArray_char_T *weight_file = allocate_for_charArray(argv[index++]);
    emxArray_char_T *motif_file  = allocate_for_charArray(argv[index++]);
    double numPWM = strtod(argv[index++], NULL);


    char arr[14][50] = {"Model prefix", "K-mer weight file", "Motif file", "Number of PWMs", "Number of iterations", 
                        "Correlation cutoff", "Regularization", "Positive-Negative PWM ratio", "K-mer fraction",
                        "Total gapped K-mer length", "Number of ungapped positions", 
                        "Average GC content", "Equivalent reverse complement K-mer", "Allow auto K-mer fraction adjustment"};
    double arr2[8] = {numPWM, numIterations, rCorr, reg, iPNRatio, kmerFrac, lSVM, kSVM};
    int arr3[3] = {backgroundGC, reverseCompl, kmerFracLimit};
    
    printf("\n=====  Following Are The Command Line Arguments Passed  =====\n");
    for(int counter=1; counter<4; counter++)
        printf("\n%s: %s", arr[counter-1], argv[optind+counter-1]);
    for(int counter=4; counter<12; counter++)
        printf("\n%s: %4.2f", arr[counter-1], arr2[counter-4]);
    for(int counter=12; counter<15; counter++)
        printf("\n%s: %s", arr[counter-1], arr3[counter-12] ? "True" : "False");
    printf("\n\n=============================================================");
    printf("\n\n");
    
#if __linux__
    openblas_set_num_threads(1);
#endif

    gkmPWM(model_file,
           weight_file,
           motif_file, 
           numPWM, 
           numIterations,
           rCorr,
           reg,
           lSVM,
           kSVM,
           backgroundGC,
           reverseCompl,
           iPNRatio, 
           kmerFrac,
           kmerFracLimit);
    
    emxFree_char_T(&model_file);    
    emxFree_char_T(&weight_file);
    emxFree_char_T(&motif_file);
    
    gkmPWM_terminate();
    return 0;
}

