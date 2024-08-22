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
#include "mapTF.h"
#include "mapTF_terminate.h"
#include "mapTF_emxAPI.h"
#include "mapTF_emxutil.h"
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
            "mapTF: A method to map TFBS motifs from gkmPWMlasso and gkmPWM to regions at base-pair resolution.\n"
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki, Gary Yang, Michael Beer\n"
            "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     mapTF [options] <sequence> <weight> <denovo motif> <lasso motif> <database> <output prefix>\n"
            "\n"
            "Arguments:\n"
            "       sequence: The set of sequences to which the motifs will be mapped (fasta format)\n"
            "         weight: file name to the two-column, tab-delimited kmer weight file \n" 
            "   denovo motif: The gkmPWM meme output file (*_denovo.meme)\n"
            "    lasso motif: The gkmPWMlasso output file (*_gkmPWMlasso.out)\n"
            "       database: The collection of position weight matrix (PWM) in meme format. \n"
            "                 MUST be the same as the one used by gkmPWMlasso\n"
            "  output prefix: output file name prefix for mapTF\n"
            "\n"
            "  Options:\n"
            "    -f <float>    Set the fraction of the total number of gapped k-mers to use with gkmPWM. \n"
            "                  This reduces the memory and runtime needed. If the total number of gapped k-mers \n"
            "                  is too high with the given combination of (l,k,KmerFrac), KmerFrac will be \n"
            "                  automatically set to a lower value to create a more workable number of gapped k-mers \n"
            "    -l <int>      total length of the gapped k-mer. It MUST be the same l from the \n"
            "                  input weight model (default: 11)\n"
            "    -k <int>      the number of ungapped positions. It MUST be the same k from the \n"
            "                  input weight model (default:  7)\n" 
            "    -L            if set, PWM probabilities for each k-mer will be saved. Setting this option speeds \n"
            "                  up mapTF; however, around 15GB of RAM is needed. Only set this if you have enough \n"
            "                  RAM. (default: false) \n"           
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
    if(argc == 1) { display_arguments(); }

    double kmerFrac = 1;
    double lSVM = 11;
    double kSVM = 7 ;
    int LS = 0;

    char * pEnd;

    int c;
    while ((c = getopt(argc, argv, "Ll:k:f:")) != -1) {
        switch (c) {
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
            case 'L':
                LS = 1;
                break;
            default:
                fprintf(stderr, "Unknown option: -%c\n", c);
                display_arguments();
        }
    }

    if (argc - optind != 6) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the .\n");
        display_arguments();
    }
    if (lSVM < kSVM) {
        printf("ERROR: k must be a positive integer smaller than l\n");
        exit(1);
    }
    
    int index = optind;
    emxArray_char_T *seq_file      = allocate_for_charArray(argv[index++]);
    emxArray_char_T *weight_file   = allocate_for_charArray(argv[index++]);
    emxArray_char_T *denovo_file   = allocate_for_charArray(argv[index++]);
    emxArray_char_T *lasso_file    = allocate_for_charArray(argv[index++]);
    emxArray_char_T *motif_file    = allocate_for_charArray(argv[index++]);
    emxArray_char_T *output_prefix = allocate_for_charArray(argv[index++]);
    
    char arr[10][30] = {"Sequence file", "K-mer weight file", "De novo motifs", "LASSO motifs", 
                       "Motif file", "Output prefix", "K-mer fraction",
                       "Total gapped k-mer length", "Number of ungapped positions", "Save sequences"};
    double arr2[4] = {kmerFrac, lSVM, kSVM, LS};
    printf("\n=====  Following Are The Command Line Arguments Passed  =====\n");
    for(int counter=1; counter<7; counter++)
        printf("\n%s: %s", arr[counter-1], argv[optind+counter-1]);
    for(int counter=7; counter<11; counter++)
        printf("\n%s: %4.2f", arr[counter-1], arr2[counter-7]);
    printf("\n\n=============================================================");
    printf("\n\n");

#if __linux__    
    openblas_set_num_threads(1);
#endif

    mapTF(seq_file,
          weight_file, 
          denovo_file,
          lasso_file,
          motif_file,
          output_prefix,
          lSVM,
          kSVM,
          kmerFrac,
          LS);
                 
    emxFree_char_T(&seq_file);
    emxFree_char_T(&weight_file);                 
    emxFree_char_T(&denovo_file);
    emxFree_char_T(&lasso_file);    
    emxFree_char_T(&motif_file);
    emxFree_char_T(&output_prefix);
    
    mapTF_terminate();
    return 0;
}

