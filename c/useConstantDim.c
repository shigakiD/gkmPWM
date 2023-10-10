/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * useConstantDim.c
 *
 * Code generation for function 'useConstantDim'
 *
 */

/* Include files */
#include "useConstantDim.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void useConstantDim(double varargin_2_data[], const int varargin_2_size[2])
{
  int i;
  int k;
  if ((varargin_2_size[1] != 0) && (varargin_2_size[1] != 1)) {
    i = varargin_2_size[1];
    for (k = 0; k <= i - 2; k++) {
      varargin_2_data[k + 1] += varargin_2_data[k];
    }
  }
}

/* End of code generation (useConstantDim.c) */
