/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * frewind.c
 *
 * Code generation for function 'frewind'
 *
 */

/* Include files */
#include "frewind.h"
#include "fileManager.h"
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void frewind(double fileID)
{
  FILE *b_NULL;
  FILE *filestar;
  bool a;
  b_NULL = NULL;
  getfilestar(fileID, &filestar, &a);
  if ((!(filestar == b_NULL)) && (fileID != 0.0) && (fileID != 1.0) &&
      (fileID != 2.0)) {
    fseek(filestar, 0, SEEK_SET);
  }
}

/* End of code generation (frewind.c) */
