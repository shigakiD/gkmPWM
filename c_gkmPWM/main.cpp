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
            "gkmPWM: A method to extract compact and interpretable " 
            "features from gkmSVM models.\n"
            "\n\n"
            "Version:   1.0\n"
            "Code:      https://github.com/shigakiD/gkmPWM/tree/main\n"
            "Author:    Dustin Shigaki\n"
            "Contact:   Report issues to the Github page\n"
            "\n\n"
            "Usage:     gkmPWM [options] <prefix> <weight file> <database> <motif number>\n"
            "\n"
            "Arguments:\n"
            "        prefix: prefix of a gkmSVM model, where the model files are either\n"
            "                *_svseq.fa and *_svalpha.out OR *.model.txt\n"
            "   weight file: file name to the corresponding ungapped k-mer weight file\n"
            "      database: file name to a database of transcription factor binding sites\n"
            "                that is in meme format. An example file is provided.\n"
            "  motif number: number of position weight matrix (PWM) to extract\n"
            "\n"
            "  Options:\n"
            "    -c <float>    correlation cutoff for considering two PWMs to be distinct (default: 0.90)\n"
            "    -r <float>    regularization (default: 0)\n"
            "    -P <float>    number of motifs to extract from the negative set used to train a gkmSVM model. \n"
            "                  when set to 0, gkmPWM will automatically infer this number. (default: 0)\n"
            "    -l <int>      total length of the gapped k-mer. It doesn't have to be the same l from gkmSVM model (default: 10)\n"
            "    -k <int>      the number of ungapped positions. It doesn't have to be the same k from gkmSVM model (default: 6)\n"
            "    -i <int>      number of training iterations. Recommend at least 10 times the number of motifs (default: 100)\n"
            "    -b            if set, background GC content equal to the average of the \n"
            "                  positive and negative set used to train gkmSVM; else, the GC content \n"
            "                  equal to the negative set (default: false) \n"
            "    -R            consider reverse-complement of k-mer to be the same (default: true) \n"
            "\n");
    exit(0);
}


int main(int argc, char* argv[]) {
    
	if(argc == 1) { display_arguments(); }

    double rCorr = 0.90;
    double reg = 0;
    double iPNRatio = 0;
    double lSVM = 10;
    double kSVM = 6 ;
    double numIterations = 100;
    int backgroundGC = 0;
    int reverseCompl = 1;
    char * pEnd;

	int c;
	while ((c = getopt(argc, argv, "bRc:r:P:l:k:i:")) != -1) {
		switch (c) {
            case 'b':
                backgroundGC = 1;
                break;
            case 'R':
                reverseCompl = 0;
                break;
            case 'c':
                rCorr = strtod(optarg, &pEnd);
                break;            
            case 'r':
                reg = strtod(optarg, &pEnd);
                break;
            case 'P':
                iPNRatio = strtod(optarg, &pEnd);
                break;
            case 'l':
                lSVM = strtod(optarg, &pEnd);
                break;
            case 'k':
                kSVM = strtod(optarg, &pEnd);
                break;
            case 'i':
                numIterations = strtod(optarg, &pEnd);
                break;
			default:
                fprintf(stderr, "Unknown option: -%c\n", c);
                display_arguments();
		}
	}

    if (argc - optind != 4) {
        fprintf(stderr, "Incorrect number of required arguments. Please read the .\n");
        display_arguments();
    }

	int index = optind;
    emxArray_char_T *model_file  = allocate_for_charArray(argv[index++]);
    emxArray_char_T *weight_file = allocate_for_charArray(argv[index++]);
    emxArray_char_T *motif_file  = allocate_for_charArray(argv[index++]);
    double numPWM = strtod(argv[index++], NULL);

    
    char arr[12][30] = {"model file", "weight file", "motif file", "number of PWM", "number of iteration", "correlation cutoff", "reg", "PN ratio", 
                        "total gapped k-mer length", "number of ungapped positions", "average GC", "reverse complement"};
    double arr2[7] = {numPWM, numIterations, rCorr, reg, iPNRatio, lSVM, kSVM};
    int arr3[2] = {backgroundGC, reverseCompl};
    
    printf("\n----Following Are The Command Line Arguments Passed----");
    for(int counter=1; counter<4; counter++)
        printf("\nargv[%d] - %s: %s", counter, arr[counter-1], argv[optind+counter-1]);
    for(int counter=4; counter<11; counter++)
        printf("\nargv[%d] - %s: %4.2f", counter, arr[counter-1], arr2[counter-4]);
    for(int counter=11; counter<13; counter++)
        printf("\nargv[%d] - %s: %s", counter, arr[counter-1], arr3[counter-11] ? "True" : "False");
    printf("\n");
    
    openblas_set_num_threads(1);
    
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
           iPNRatio);
    
    emxFree_char_T(&model_file);    
    emxFree_char_T(&weight_file);
    emxFree_char_T(&motif_file);
    
    gkmPWM_terminate();
    return 0;
}

// End of code generation (main.cpp)
