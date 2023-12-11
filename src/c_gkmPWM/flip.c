/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * flip.c
 *
 * Code generation for function 'flip'
 *
 */

/* Include files */
#include "flip.h"
#include "gkmPWM_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void flip(emxArray_real_T *x)
{
  double tmp;
  double *x_data;
  int b_i;
  int dim;
  int i;
  int i1;
  int j;
  int k;
  int nd2;
  int npages;
  int offset;
  int pagelen;
  int tmp_tmp;
  int vstride;
  x_data = x->data;
  dim = 1;
  if (x->size[0] != 1) {
    dim = 0;
  }
  if ((x->size[0] != 0) && (x->size[1] != 0) && (x->size[dim] > 1)) {
    vstride = 1;
    for (k = 0; k < dim; k++) {
      vstride *= x->size[0];
    }
    pagelen = vstride * x->size[dim];
    npages = 1;
    i = dim + 2;
    for (k = i; k < 3; k++) {
      npages *= x->size[1];
    }
    i = x->size[dim] - 1;
    nd2 = x->size[dim] >> 1;
    dim = npages - 1;
    for (j = 0; j <= dim; j++) {
      npages = vstride - 1;
      for (b_i = 0; b_i <= npages; b_i++) {
        offset = j * pagelen + b_i;
        for (k = 0; k < nd2; k++) {
          tmp_tmp = offset + k * vstride;
          tmp = x_data[tmp_tmp];
          i1 = offset + (i - k) * vstride;
          x_data[tmp_tmp] = x_data[i1];
          x_data[i1] = tmp;
        }
      }
    }
  }
}

/* End of code generation (flip.c) */
