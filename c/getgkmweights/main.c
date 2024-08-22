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
#include "getgkmweights.h"
#include "getgkmweights_terminate.h"
#include "getgkmweights_emxAPI.h"
#include "getgkmweights_emxutil.h"
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
            "getgkmweights: A method to generate the gapped kmers weights for a (l,k) pair\n"
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki, Gary Yang, Michael Beer\n"
            "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     getgkmweights [options] <prefix> <l_svm> <k_svm>\n"
            "\n"
            "Arguments:\n"
            "                     prefix: prefix of a gkmSVM model, where the model files are either\n"
            "                             FILEHEADER_svseq.fa and FILEHEADER_svalpha.out OR FILEHEADER.model.txt\n"
            "                      l_svm: The full length of the gapped k-mer.  This DOES NOT need to be the same as the \n"
            "                             l in the gkmSVM model. \n"
            "                      k_svm: The number of ungapped positions of the gapped k-mer. This DOES NOT need to be \n"
            "                             the same as the k in the gkmSVM model. \n"
            "\n"
            "  Options:\n"
            "    -f <float>    Set the fraction of the total number of gapped k-mers to use with gkmPWM. \n"
            "                  This reduces the memory and runtime needed. If the total number of gapped k-mers \n"
            "                  is too high with the given combination of (l,k,KmerFrac), KmerFrac will be \n"
            "                  automatically set to a lower value to create a more workable number of gapped k-mers. \n"
            "                  (default: 1) \n"
            "    -R            if set, consider reverse-complements of gapped k-mers to be distinct features \n"
            "    -L            if set, WILL NOT automatically lower the fraction of gapped-kmer used. When set, \n"
            "                  the program will require more memory and need longer runtime \n"
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
    if(argc == 1) { display_arguments(); }

    double kmerFrac = 1;
    int reverseCompl = 1;
    int kmerFracLimit = 1;
    char * pEnd;

    int c;
    while ((c = getopt(argc, argv, "LRf:")) != -1) {
        switch (c) {
            case 'R':
                reverseCompl = 0;
                break;      
             case 'L':
                kmerFracLimit = 0;
                break;
             case 'f':
                kmerFrac = strtod(optarg, &pEnd);
                if (kmerFrac <= 0 || kmerFrac > 1) {
                    printf("ERROR: Fraction K-mer used be a positive value in (0, 1].\n");
                    exit(1);
                }
                break;
            default:
                fprintf(stderr, "Unknown option: -%c\n", c);
                display_arguments();
        }
    }

    if (argc - optind != 3) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the documentation.\n");
        display_arguments();
    }
    int index = optind;
    emxArray_char_T *model_file  = allocate_for_charArray(argv[index++]);
    double lSVM = strtod(argv[index++], &pEnd);
    double kSVM = strtod(argv[index++], &pEnd);
    
    if (lSVM < kSVM || lSVM <= 0 || kSVM <= 0) {
        printf("ERROR: k must be a positive integer smaller than l\n");
        exit(1);
    }
    

    char arr[6][50] = {"Model prefix", "Total gapped K-mer length", "Number of ungapped positions", 
                       "K-mer fraction", "Equivalent reverse complement K-mers", 
                       "Allow auto K-mer fraction adjustment"};
    double arr2[3] = {kmerFrac};
    int arr3[2] = {reverseCompl, kmerFracLimit};
    
    printf("\n=====  Following Are The Command Line Arguments Passed  =====\n");
    for(int counter=1; counter<4; counter++)
        printf("\n%s: %s", arr[counter-1], argv[optind+counter-1]);
    for(int counter=4; counter<5; counter++)
        printf("\n%s: %4.2f", arr[counter-1], arr2[counter-4]);
    for(int counter=5; counter<7; counter++)
        printf("\n%s: %s", arr[counter-1], arr3[counter-5] ? "True" : "False");
    printf("\n\n=============================================================");
    printf("\n\n");
    
#if __linux__
    openblas_set_num_threads(1);
#endif

    getgkmweights(model_file, lSVM, kSVM, reverseCompl, kmerFrac, kmerFracLimit);
    
    emxFree_char_T(&model_file);    
    
    getgkmweights_terminate();
    return 0;
}

