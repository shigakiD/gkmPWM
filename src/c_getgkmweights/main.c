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
#include "lapacke.h"
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
            "Usage:     gkmPWM [options] <prefix> \n"
            "\n"
            "Arguments:\n"
            "                     prefix: prefix of a gkmSVM model, where the model files are either\n"
            "                             FILEHEADER_svseq.fa and FILEHEADER_svalpha.out OR FILEHEADER.model.txt\n"
            "\n"
            "  Options:\n"
	    "    -f <float>    Set the fraction of the total number of gapped k-mers to use with gkmPWM. \n"
            "                  This reduces the memory and runtime needed. If the total number of gapped k-mers \n"
            "                  is too high with the given combination of (l,k,KmerFrac), KmerFrac will be \n"
	    "                  automatically set to a lower value to create a more workable number of gapped k-mers \n"
            "    -l <int>      The full length of the gapped k-mer.  This DOES NOT need to be the same as the \n"
	    "                  l in the gkmSVM model (default: 10) \n"
            "    -k <int>      The number of ungapped positions of the gapped k-mer. This DOES NOT need to be \n"
	    "                  the same as the k in the gkmSVM model (default: 6)\n"
            "    -R            if set, consider reverse-complements of gapped k-mers to be distinct \n"
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
	if(argc == 1) { display_arguments(); }

    double lSVM = 10;
    double kSVM = 6 ;
    double kmerFrac = 1;
    int reverseCompl = 1;
    char * pEnd;

	int c;
	while ((c = getopt(argc, argv, "Rl:k:f:")) != -1) {
		switch (c) {
            case 'R':
                reverseCompl = 0;
                break;      
            case 'l':
                lSVM = strtod(optarg, &pEnd);
                break;
            case 'k':
                kSVM = strtod(optarg, &pEnd);
                break;
            case 'f':
                kmerFrac = strtod(optarg, &pEnd);
                break;
			default:
                fprintf(stderr, "Unknown option: -%c\n", c);
                display_arguments();
		}
	}

    if (argc - optind != 1) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the documentation.\n");
        display_arguments();
    }

	int index = optind;
    emxArray_char_T *model_file  = allocate_for_charArray(argv[index++]);

    
    char arr[5][30] = {"model file",  "kmer fraction",
                        "total gapped k-mer length", "number of ungapped positions", "reverse complement"};
    double arr2[3] = {kmerFrac, lSVM, kSVM};
    int arr3[1] = {reverseCompl};
    
    printf("\n----Following Are The Command Line Arguments Passed----");
    for(int counter=1; counter<2; counter++)
        printf("\nargv[%d] - %s: %s", counter, arr[counter-1], argv[optind+counter-1]);
    for(int counter=2; counter<5; counter++)
        printf("\nargv[%d] - %s: %4.2f", counter, arr[counter-1], arr2[counter-2]);
    for(int counter=5; counter<6; counter++)
        printf("\nargv[%d] - %s: %s", counter, arr[counter-1], arr3[counter-5] ? "True" : "False");
    printf("\n");
    
    openblas_set_num_threads(1);

    getgkmweights(model_file, lSVM, kSVM, reverseCompl, kmerFrac);
    
    emxFree_char_T(&model_file);    
    
    getgkmweights_terminate();
    return 0;
}

// End of code generation (main.cpp)
