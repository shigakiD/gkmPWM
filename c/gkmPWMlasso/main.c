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
#include "gkmPWMlasso.h"
#include "gkmPWMlasso_terminate.h"
#include "gkmPWMlasso_emxAPI.h"
#include "gkmPWMlasso_emxutil.h"
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
            "gkmPWMlasso: A method to extract compact predictive motifs from sequence-based models of regulatory elements using a database of PWMs (meme format).\n" 
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki, Gary Yang, Michael Beer\n"
            "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     gkmPWMlasso [options] <prefix> <database> <motif number>\n"
            "\n"
            "Arguments:\n"
            "       prefix: prefix of a gkmSVM model, where the model files are either\n"
            "               FILEHEADER_svseq.fa and FILEHEADER_svalpha.out OR FILEHEADER.model.txt\n"
            "     database: file name to a database of transcription factor binding sites\n"
            "               that is in meme format. An example file (combined_db_v4.meme) is provided.\n"
            " motif number: The number of motifs to learn.  If 0 is specified, gkmPWMlasso will \n"
            "               try to find the optimal number of motifs. \n"
            "\n"
            "  Options:\n"
            "    -m <int>      PWMs with lengths shorter than this will be filtered out (default: 10)\n"
            "    -i <float>    PWMs with an average information per position less than this will be filtered out (default: 0.5)\n"
            "    -c <float>    PWMs with a pearson correlation greater than this will be clustered \n"
            "                  together to prevent linear dependence and redundant features (default: 0.83)\n"
            "    -f <float>    Set the fraction of the total number of gapped k-mers to use with gkmPWM. \n"
            "                  This reduces the memory and runtime needed. If the total number of gapped k-mers \n"
            "                  is too high with the given combination of (l,k,KmerFrac), KmerFrac WILL be \n"
            "                  automatically set to a lower value to create a more workable number of gapped k-mers \n"
            "    -l <int>      The full length of the gapped k-mer. This DOES NOT need to be the same as the \n"
            "                  l in the gkmSVM model (default: 11) \n"
            "    -k <int>      The number of ungapped positions of the gapped k-mer. This DOES NOT need to be \n"
            "                  the same as the k in the gkmSVM model (default: 7) \n"
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

    double minLength = 10;
    double minInfo = 0.5;
    double corrCut = 0.83;
    double kmerFrac = 1;
    double lSVM = 11;
    double kSVM = 7 ;
    int backgroundGC = 0;
    int reverseCompl = 1;
    int kmerFracLimit = 1;
    char * pEnd;

    int c;
    while ((c = getopt(argc, argv, "BRLd:m:i:c:l:k:f:")) != -1) {
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
            case 'm':
                minLength = strtod(optarg, &pEnd);
                if (minLength <= 0) {
                    printf("ERROR: Minimum length must be a positive integer.\n");
                    exit(1);
                }
                break;
            case 'i':
                minInfo = strtod(optarg, &pEnd);
                if (minInfo <= 0) {
                    printf("ERROR: Minimum information content must be a positive value.\n");
                    exit(1);
                }
                break;
            case 'c':
                corrCut = strtod(optarg, &pEnd);
                if (corrCut <= 0) {
                    printf("ERROR: Correlation cutoff must be a positive value.\n");
                    exit(1);
                }
                break;
            case 'f':
                kmerFrac = strtod(optarg, &pEnd);
                if (kmerFrac <= 0 || kmerFrac > 1) {
                    printf("ERROR: Fraction K-mer used be a positive value in (0, 1].\n");
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
            default:
                fprintf(stderr, "Unknown option: -%c\n", c);
                display_arguments();
        }
    }

    if (argc - optind != 3) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the .\n");
        display_arguments();
    }
    if (lSVM < kSVM) {
        printf("ERROR: k must be a positive integer smaller than l\n");
        exit(1);
    }

    int index = optind;
    emxArray_char_T *model_file = allocate_for_charArray(argv[index++]);
    emxArray_char_T *motif_file = allocate_for_charArray(argv[index++]);
    double numPWM = strtod(argv[index++], &pEnd);

    char arr[12][50] = {"Model prefix", "Motif file", "Number of PWMs", "Minimum PWM length", 
                "Minimum PWM information content", "Correlation cutoff", "K-mer fraction",
                "Total gapped K-mer length", "Number of ungapped positions", 
                "Average GC content", "Equivalent reverse complement K-mer", 
                "Allow auto K-mer fraction adjustment"};
    double arr2[6] = {minLength, minInfo, corrCut, kmerFrac, lSVM, kSVM};
    int arr3[3] = {backgroundGC, reverseCompl, kmerFracLimit};
    printf("\n=====  Following Are The Command Line Arguments Passed  =====\n");
    for(int counter=1; counter<4; counter++)
        printf("\n%s: %s", arr[counter-1], argv[optind+counter-1]);
    for(int counter=4; counter<10; counter++)
        printf("\n%s: %4.2f", arr[counter-1], arr2[counter-4]);
    for(int counter=10; counter<13; counter++)
        printf("\n%s: %s", arr[counter-1], arr3[counter-10] ? "True" : "False");
    printf("\n\n=============================================================");
    printf("\n\n");

#if __linux__
    openblas_set_num_threads(1);
#endif
    
    gkmPWMlasso(model_file,
                motif_file, numPWM,
                minLength,
                minInfo,
                corrCut,
                lSVM,
                kSVM,
                backgroundGC,
                reverseCompl,
                kmerFrac, 
                kmerFracLimit);
    
    emxFree_char_T(&model_file);
    emxFree_char_T(&motif_file);
    
    gkmPWMlasso_terminate();
    return 0;
}

