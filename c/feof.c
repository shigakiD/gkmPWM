/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * feof.c
 *
 * Code generation for function 'feof'
 *
 */

/* Include files */
#include "feof.h"
#include "fileManager.h"
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
double b_feof(double fileID)
{
  FILE *b_NULL;
  FILE *filestar;
  int b_st;
  double st;
  bool a;
  b_NULL = NULL;
  getfilestar(fileID, &filestar, &a);
  if (filestar == b_NULL) {
    st = 0.0;
  } else {
    b_st = feof(filestar);
    st = ((int)b_st != 0);
  }
  return st;
}

/* End of code generation (feof.c) */
