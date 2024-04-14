/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fseek.c
 *
 * Code generation for function 'fseek'
 *
 */

/* Include files */
#include "fseek.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_fseek(double fileID, double offset)
{
  FILE *filestar;
  int wherefrom;
  bool a;
  wherefrom = SEEK_SET;
  if ((!rtIsInf(offset)) && (!rtIsNaN(offset)) && (floor(offset) == offset)) {
    getfilestar(fileID, &filestar, &a);
    if ((!(fileID != 0.0)) || (!(fileID != 1.0)) || (!(fileID != 2.0))) {
      filestar = NULL;
    }
    if (!(filestar == NULL)) {
      fseek(filestar, (long int)offset, wherefrom);
    }
  }
}

/* End of code generation (fseek.c) */
