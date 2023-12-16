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
            "gkmPWMlasso: A method to extract compact predictive motifs from sequence-based models of regulatory elements using a database of PWMs (meme format)." 
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki, Gary Yang, Michael Beer\n"
	    "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     gkmPWMlasso [options] <prefix> <database> <motif number>"
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
	    "                  together to prevent linear dependence and redundant features (default: 0.86)\n"
	    "    -f <float>    Set the fraction of the total number of gapped k-mers to use with gkmPWM. \n"
            "                  This reduces the memory and runtime needed. If the total number of gapped k-mers \n"
            "                  is too high with the given combination of (l,k,KmerFrac), KmerFrac will be \n"
	    "                  automatically set to a lower value to create a more workable number of gapped k-mers \n"
            "    -l <int>      The full length of the gapped k-mer. This DOES NOT need to be the same as the \n"
	    "                  l in the gkmSVM model (default: 10) \n"
            "    -k <int>      The number of ungapped positions of the gapped k-mer. This DOES NOT need to be \n"
	    "                  the same as the k in the gkmSVM model (default: 6) \n"
            "    -b            if set, use both the positive and negative set to get the background \n"
            "                  distribution of gapped k-mers. This is best used when the model is trained \n"
            "                  with both the positive and negative sets containing regulatory elements.\n"
            "    -R            if set, consider reverse-complements of gapped k-mers to be distinct \n"
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
	if(argc == 1) { display_arguments(); }

	double minLength = 10;
	double minInfo = 0.5;
	double corrCut = 0.86;
	double kmerFrac = 1;
	double lSVM = 10;
	double kSVM = 6 ;
	int backgroundGC = 0;
	int reverseCompl = 1;
	char * pEnd;

	int c;
	while ((c = getopt(argc, argv, "bRd:m:i:c:l:k:f:")) != -1) {
		switch (c) {
            case 'b':
                backgroundGC = 1;
                break;
            case 'R':
                reverseCompl = 0;
                break;
            case 'm':
                minLength = strtod(optarg, &pEnd);
                break;
            case 'i':
                minInfo = strtod(optarg, &pEnd);
                break;
            case 'c':
                corrCut = strtod(optarg, &pEnd);
                break;
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

	if (argc - optind != 3) {
		fprintf(stderr, "Incorrect number of required arguments. Please read the .\n");
		display_arguments();
	}

	int index = optind;
	emxArray_char_T *model_file = allocate_for_charArray(argv[index++]);
	emxArray_char_T *motif_file = allocate_for_charArray(argv[index++]);
	double numPWM = strtod(argv[index++], &pEnd);

	char arr[11][30] = {"model file", "motif file", "number of PWM", "min PWM length", 
				"min PWM info", "correlation cutoff", "kmer fraction",
				"total gapped k-mer length", "number of ungapped positions", 
				"average GC", "reverse complement"};
	double arr2[6] = {minLength, minInfo, corrCut, kmerFrac, lSVM, kSVM};
	int arr3[2] = {backgroundGC, reverseCompl};
	printf("\n----Following Are The Command Line Arguments Passed----");
	for(int counter=1; counter<4; counter++)
		printf("\nargv[%d] - %s: %s", counter, arr[counter-1], argv[optind+counter-1]);
	for(int counter=4; counter<10; counter++)
		printf("\nargv[%d] - %s: %4.2f", counter, arr[counter-1], arr2[counter-4]);
	for(int counter=10; counter<12; counter++)
		printf("\nargv[%d] - %s: %s", counter, arr[counter-1], arr3[counter-10] ? "True" : "False");
	printf("\n");
    
	openblas_set_num_threads(1);
    
	gkmPWMlasso(model_file,
			motif_file, 
			minLength,
			minInfo,
			corrCut,
			lSVM,
			kSVM,
			backgroundGC,
			reverseCompl,
			numPWM, 
			kmerFrac);
    
	emxFree_char_T(&model_file);
	emxFree_char_T(&motif_file);
    
	gkmPWMlasso_terminate();
	return 0;
}

// End of code generation (main.c)
