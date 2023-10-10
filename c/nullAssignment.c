/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nullAssignment.c
 *
 * Code generation for function 'nullAssignment'
 *
 */

/* Include files */
#include "nullAssignment.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_nullAssignment(emxArray_real_T *x, const emxArray_int32_T *idx)
{
  emxArray_boolean_T *b;
  double *x_data;
  const int *idx_data;
  int i;
  int k;
  int k0;
  int nxin;
  int nxout;
  bool *b_data;
  idx_data = idx->data;
  x_data = x->data;
  emxInit_boolean_T(&b, 2);
  nxin = x->size[0];
  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = x->size[0];
  emxEnsureCapacity_boolean_T(b, i);
  b_data = b->data;
  k0 = x->size[0];
  for (i = 0; i < k0; i++) {
    b_data[i] = false;
  }
  i = idx->size[0];
  for (k = 0; k < i; k++) {
    b_data[idx_data[k] - 1] = true;
  }
  k0 = 0;
  i = b->size[1];
  for (k = 0; k < i; k++) {
    k0 += b_data[k];
  }
  nxout = x->size[0] - k0;
  k0 = -1;
  for (k = 0; k < nxin; k++) {
    if ((k + 1 > b->size[1]) || (!b_data[k])) {
      k0++;
      x_data[k0] = x_data[k];
    }
  }
  emxFree_boolean_T(&b);
  i = x->size[0];
  if (1 > nxout) {
    x->size[0] = 0;
  } else {
    x->size[0] = nxout;
  }
  emxEnsureCapacity_real_T(x, i);
}

/*
 *
 */
void c_nullAssignment(emxArray_real_T *x, int idx)
{
  double *x_data;
  int i;
  int j;
  int ncols;
  int ncolx;
  int nrowx;
  x_data = x->data;
  nrowx = x->size[0];
  ncolx = x->size[1] - 2;
  ncols = x->size[1] - 1;
  for (j = idx; j <= ncols; j++) {
    for (i = 0; i < nrowx; i++) {
      x_data[i + x->size[0] * (j - 1)] = x_data[i + x->size[0] * j];
    }
  }
  if (1 > ncols) {
    nrowx = -1;
  } else {
    nrowx = ncolx;
  }
  ncolx = x->size[0] - 1;
  ncols = x->size[0];
  for (j = 0; j <= nrowx; j++) {
    for (i = 0; i < ncols; i++) {
      x_data[i + (ncolx + 1) * j] = x_data[i + x->size[0] * j];
    }
  }
  j = x->size[0] * x->size[1];
  x->size[0] = ncolx + 1;
  x->size[1] = nrowx + 1;
  emxEnsureCapacity_real_T(x, j);
}

/*
 *
 */
void nullAssignment(emxArray_real_T *x, int idx)
{
  double *x_data;
  int i;
  int j;
  int ncolx;
  int nrows;
  int nrowx;
  x_data = x->data;
  nrowx = x->size[0] - 2;
  ncolx = x->size[1];
  nrows = x->size[0] - 1;
  for (j = 0; j < ncolx; j++) {
    for (i = idx; i <= nrows; i++) {
      x_data[(i + x->size[0] * j) - 1] = x_data[i + x->size[0] * j];
    }
  }
  if (1 > nrows) {
    nrowx = 0;
  } else {
    nrowx++;
  }
  ncolx = x->size[1] - 1;
  for (nrows = 0; nrows <= ncolx; nrows++) {
    for (j = 0; j < nrowx; j++) {
      x_data[j + nrowx * nrows] = x_data[j + x->size[0] * nrows];
    }
  }
  nrows = x->size[0] * x->size[1];
  x->size[0] = nrowx;
  x->size[1] = ncolx + 1;
  emxEnsureCapacity_real_T(x, nrows);
}

/* End of code generation (nullAssignment.c) */
