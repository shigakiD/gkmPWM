/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fliplr.c
 *
 * Code generation for function 'fliplr'
 *
 */

/* Include files */
#include "fliplr.h"
#include "gkmPWMlasso3_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void fliplr(emxArray_real_T *x)
{
  double xtmp;
  double *x_data;
  int b_j1;
  int j2;
  int n;
  int nd2;
  x_data = x->data;
  n = x->size[1] - 1;
  nd2 = x->size[1] >> 1;
  for (b_j1 = 0; b_j1 < nd2; b_j1++) {
    j2 = n - b_j1;
    xtmp = x_data[b_j1];
    x_data[b_j1] = x_data[j2];
    x_data[j2] = xtmp;
  }
}

/* End of code generation (fliplr.c) */
