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
            "mapTF: A method to map motifs from gkmPWMlasso and gkmPWM to regions at base-pair resolution.\n"
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki, Gary Yang, Michael Beer\n"
	    "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     mapTF [options] <sequence> <weight> <denovo motif> <lasso motif> <database> <output prefix>"
            "\n"
            "Arguments:\n"
            "       sequence: sequence file (fasta format)\n"
            "         weight: file name to the ungapped k-mer weight file\n" 
            "   denovo motif: file name to the motifs found by gkmPWM (*_denovo.meme)\n"
            "    lasso motif: file name to the output generated by gkmPWMlasso (*_gkmPWMlasso.out)\n"
            "       database: file name to the database of transcription factor binding sites\n"
            "                 (meme format). MUST be the same as the one used by gkmPWMlasso\n"
            "  output prefix: output file name prefix for mapTF\n"
            "\n"
            "  Options:\n"
	    "    -f <float>    a fraction of total gapped k-mer for PWM extraction; the number must be in [0,1];\n"
            "                  setting to a smaller value, like 0.2, will significantly speed up the program.\n"
            "                  If the total number of gapped k-mer is too large, this value will be set \n"
	    "                  automatically to a smaller number (default: 1)\n"
            "    -l <int>      total length of the gapped k-mer. It MUST be the same l from the \n"
	    "                  input weight model (default: 11)\n"
            "    -k <int>      the number of ungapped positions. It MUST be the same k from the \n"
	    "                  input weight model (default:  7)\n" 
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
	if(argc == 1) { display_arguments(); }

    double kmerFrac = 1;
    double lSVM = 11;
    double kSVM = 7 ;

	char * pEnd;

	int c;
	while ((c = getopt(argc, argv, "l:k:f:")) != -1) {
		switch (c) {
            case 'f':
                kmerFrac = strtod(optarg, &pEnd);
                break;
            case 'l':
                lSVM = strtod(optarg, &pEnd);
                break;
            case 'k':
                kSVM = strtod(optarg, &pEnd);
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

	int index = optind;
    emxArray_char_T *seq_file      = allocate_for_charArray(argv[index++]);
    emxArray_char_T *weight_file   = allocate_for_charArray(argv[index++]);
    emxArray_char_T *denovo_file   = allocate_for_charArray(argv[index++]);
    emxArray_char_T *lasso_file    = allocate_for_charArray(argv[index++]);
    emxArray_char_T *motif_file    = allocate_for_charArray(argv[index++]);
    emxArray_char_T *output_prefix = allocate_for_charArray(argv[index++]);
    
    char arr[9][30] = {"sequence file", "weight file", "de novo motifs", "lasso motifs", 
                        "motif database", "output prefix", "kmer fraction",
                        "total gapped k-mer length", "number of ungapped positions"};
    double arr2[3] = {kmerFrac, lSVM, kSVM};
    printf("\n----Following Are The Command Line Arguments Passed----");
    for(int counter=1; counter<7; counter++)
        printf("\nargv[%d] - %s: %s", counter, arr[counter-1], argv[optind+counter-1]);
    for(int counter=7; counter<10; counter++)
        printf("\nargv[%d] - %s: %4.2f", counter, arr[counter-1], arr2[counter-7]);

    printf("\n");
    
    openblas_set_num_threads(1);
    
    mapTF(seq_file,
                 weight_file, 
                 denovo_file,
                 lasso_file,
                 motif_file,
                 output_prefix,
                 lSVM,
                 kSVM,
                 kmerFrac);
                 
    emxFree_char_T(&seq_file);
    emxFree_char_T(&weight_file);                 
    emxFree_char_T(&denovo_file);
    emxFree_char_T(&lasso_file);    
    emxFree_char_T(&motif_file);
    emxFree_char_T(&output_prefix);
    
    mapTF_terminate();
    return 0;
}

// End of code generation (main.cpp)