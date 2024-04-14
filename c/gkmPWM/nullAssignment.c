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
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_nullAssignment(emxArray_real_T *x, int idx)
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
void c_nullAssignment(double x_data[], int x_size[2], const int idx[2])
{
  int k;
  int k0;
  int n;
  bool b_data[4];
  bool b;
  for (n = 0; n < 4; n++) {
    b_data[n] = false;
  }
  for (k = 0; k < 2; k++) {
    b_data[idx[k] - 1] = true;
  }
  n = 0;
  k0 = -1;
  for (k = 0; k < 4; k++) {
    b = b_data[k];
    n += b;
    if (!b) {
      k0++;
      x_data[k0] = x_data[k];
    }
  }
  if (1 > 4 - n) {
    x_size[1] = 0;
  } else {
    x_size[1] = 4 - n;
  }
}

/*
 *
 */
void d_nullAssignment(double x_data[], int x_size[2], const int idx[2])
{
  int i;
  int j;
  int k;
  int n;
  bool b_data[4];
  bool b;
  for (n = 0; n < 4; n++) {
    b_data[n] = false;
  }
  for (k = 0; k < 2; k++) {
    b_data[idx[k] - 1] = true;
  }
  n = 0;
  j = 0;
  for (k = 0; k < 4; k++) {
    b = b_data[k];
    n += b;
    if (!b) {
      for (i = 0; i < 4; i++) {
        x_data[i + 4 * j] = x_data[i + 4 * k];
      }
      j++;
    }
  }
  x_size[0] = 4;
  if (1 > 4 - n) {
    j = -1;
  } else {
    j = 3 - n;
  }
  x_size[1] = j + 1;
}

/*
 *
 */
void e_nullAssignment(double x_data[], int x_size[2], const int idx[2])
{
  int b_i;
  int i;
  int j;
  int k;
  int n;
  int ncolx;
  bool b_data[4];
  bool b;
  ncolx = x_size[1];
  for (i = 0; i < 4; i++) {
    b_data[i] = false;
  }
  for (k = 0; k < 2; k++) {
    b_data[idx[k] - 1] = true;
  }
  n = 0;
  b_i = 0;
  for (k = 0; k < 4; k++) {
    b = b_data[k];
    n += b;
    if (!b) {
      for (j = 0; j < ncolx; j++) {
        i = x_size[0] * j;
        x_data[b_i + i] = x_data[k + i];
      }
      b_i++;
    }
  }
  if (1 > 4 - n) {
    ncolx = 0;
  } else {
    ncolx = 4 - n;
  }
  n = x_size[1] - 1;
  for (i = 0; i <= n; i++) {
    for (b_i = 0; b_i < ncolx; b_i++) {
      x_data[b_i + ncolx * i] = x_data[b_i + x_size[0] * i];
    }
  }
  x_size[0] = ncolx;
  x_size[1] = n + 1;
}

/*
 *
 */
void nullAssignment(emxArray_real_T *x, const emxArray_int32_T *idx)
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
  nxin = x->size[1];
  i = b->size[0] * b->size[1];
  b->size[0] = 1;
  b->size[1] = x->size[1];
  emxEnsureCapacity_boolean_T(b, i);
  b_data = b->data;
  k0 = x->size[1];
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
  nxout = x->size[1] - k0;
  k0 = -1;
  for (k = 0; k < nxin; k++) {
    if ((k + 1 > b->size[1]) || (!b_data[k])) {
      k0++;
      x_data[k0] = x_data[k];
    }
  }
  emxFree_boolean_T(&b);
  i = x->size[0] * x->size[1];
  if (1 > nxout) {
    x->size[1] = 0;
  } else {
    x->size[1] = nxout;
  }
  emxEnsureCapacity_real_T(x, i);
}

/* End of code generation (nullAssignment.c) */
