//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// main.cpp
//
// Code generation for function 'main'
//

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

// Include files
#include "main.h"
#include "gkmPWMlasso3.h"
#include "gkmPWMlasso3_terminate.h"
#include "gkmPWMlasso3_emxAPI.h"
#include "gkmPWMlasso3_emxutil.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>  // for strtol

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
  
int main(int argc, char* argv[]) {
    int counter;
    printf("Program Name Is: %s\n",argv[0]);
    if(argc!=10)
        printf("Missing Argument. Please specify all 9 parameters.\n");
    if(argc==10) {
        printf("\nNumber Of Arguments Passed: %d",argc);
        printf("\n----Following Are The Command Line Arguments Passed----");
        for(counter=0;counter<argc;counter++)
            printf("\nargv[%d]: %s",counter,argv[counter]);
	    char * pEnd;
	    printf("\n");
	    emxArray_char_T *model_file = allocate_for_charArray(argv[1]);
	    emxArray_char_T *motif_file = allocate_for_charArray(argv[2]);
	
	    gkmPWMlasso3(model_file,
	                 motif_file, 
	                 strtod(argv[3], &pEnd), 		         
	                 strtod(argv[4], &pEnd),
		             strtod(argv[5], &pEnd), 
		             strtod(argv[6], &pEnd), 
		             strtod(argv[7], &pEnd), 
		             strtod(argv[8], &pEnd),
		             strtod(argv[9], &pEnd));
	    emxFree_char_T(&model_file);
	    emxFree_char_T(&motif_file);
    }
    gkmPWMlasso3_terminate();
    return 0;
}

// End of code generation (main.cpp)
