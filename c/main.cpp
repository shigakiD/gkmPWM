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
#include "gkmPWMlasso4.h"
#include "gkmPWMlasso4_terminate.h"
#include "gkmPWMlasso4_emxAPI.h"
#include "gkmPWMlasso4_emxutil.h"
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
            "gkmPWMlasso: A method to extract compact and interpretable " 
            "features from gkmSVM models.\n"
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Contact:   dshigak1@jhmi.edu\n"
            "\n\n"
            "Usage:     gkmPWMlasso [options] <prefix> <database>"
            "\n"
            "Arguments:\n"
            "    prefix: prefix of a gkmSVM model, where the model files are either\n"
            "            *_svseq.fa and *_svalpha.out or\n"
            "            *.model.txt\n"
            "  database: file name to a database of transcription factor binding sites\n"
            "            that is in meme format. An example file is provided.\n"
            "\n"
            "  Options:\n"
            "    -d <int>      number of position weight matrix (PWM) to extract from gkmSVM model, \n"
            "                  if set to 0, gkmPWMlasso will automatically determine the number of PWMs \n"
            "                  that explains 90 percent of explainable variance"
            "    -m <int>      minimum accepted length of a database PWM (default: 10)\n"
            "    -i <float>    minimum accepted information content of a database PWM (default: 0.5)\n"
            "    -c <float>    correlation cutoff for considering two PWMs to be distinct (default: 0.86)\n"
            "    -l <int>      length of ungapped k-mer used when training gkmSVM (default: 10)\n"
            "    -k <int>      length of   gapped k-mer used when training gkmSVM (default: 6)\n"
            "    -b            if set, background GC content equal to the average of the \n"
            "                  positive and negative set used to train gkmSVM; else, the GC content \n"
            "                  equal to the negative set\n"
            "    -R            if set, reverse-complement of k-mer is not considered as the same feature \n"
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
	if(argc == 1) { display_arguments(); }

    double numPWM = 0;
    double minLength = 10;
    double minInfo = 0.5;
    double corrCut = 0.86;
    double lSVM = 10;
    double kSVM = 6 ;
    int backgroundGC = 0;
    int reverseCompl = 1;
	char * pEnd;

	int c;
	while ((c = getopt(argc, argv, "bRd:m:i:c:l:k:")) != -1) {
		switch (c) {
            case 'b':
                backgroundGC = 1;
                break;
            case 'R':
                reverseCompl = 1;
                break;
            case 'd':
                numPWM = strtod(optarg, &pEnd);
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

    if (argc - optind != 2) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the .\n");
        display_arguments();
    }

	int index = optind;
    emxArray_char_T *model_file = allocate_for_charArray(argv[index++]);
    emxArray_char_T *motif_file = allocate_for_charArray(argv[index++]);

    char arr[10][30] = {"model file", "motif file", "number of PWM", "min PWM length", "min PWM info", "correlation cutoff", 
                        "ungapped k-mer length", "gapped k-mer length", "average GC", "reverse complement"};
    double arr2[6] = {numPWM, minLength, minInfo, corrCut, lSVM, kSVM};
    int arr3[2] = {backgroundGC, reverseCompl};
    printf("\n----Following Are The Command Line Arguments Passed----");
    for(int counter=1; counter<3; counter++)
        printf("\nargv[%d] - %s: %s", counter, arr[counter-1], argv[optind+counter-1]);
    for(int counter=2; counter<8; counter++)
        printf("\nargv[%d] - %s: %4.2f", counter+1, arr[counter], arr2[counter-2]);
    for(int counter=8; counter<10; counter++)
        printf("\nargv[%d] - %s: %s", counter+1, arr[counter], arr3[counter-8] ? "True" : "False");
    printf("\n");
    
    openblas_set_num_threads(1);
    
    gkmPWMlasso4(model_file,
                 motif_file, 
                 minLength,
                 minInfo,
                 corrCut,
                 lSVM,
                 kSVM,
                 backgroundGC,
                 reverseCompl,
                 numPWM);
    
    emxFree_char_T(&model_file);
    emxFree_char_T(&motif_file);
    
    gkmPWMlasso4_terminate();
    return 0;
}

// End of code generation (main.cpp)
