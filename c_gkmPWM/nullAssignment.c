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
void b_nullAssignment(double x_data[], int x_size[2], const int idx[2])
{
  int j;
  int k;
  int n;
  bool b_data[4];
  bool b;
  b_data[0] = false;
  b_data[1] = false;
  b_data[2] = false;
  b_data[3] = false;
  b_data[idx[0] - 1] = true;
  b_data[idx[1] - 1] = true;
  n = 0;
  j = 0;
  for (k = 0; k < 4; k++) {
    b = b_data[k];
    n += b;
    if (!b) {
      x_data[4 * j] = x_data[4 * k];
      x_data[4 * j + 1] = x_data[4 * k + 1];
      x_data[4 * j + 2] = x_data[4 * k + 2];
      x_data[4 * j + 3] = x_data[4 * k + 3];
      j++;
    }
  }
  x_size[0] = 4;
  if (1 > 4 - n) {
    n = -1;
  } else {
    n = 3 - n;
  }
  x_size[1] = n + 1;
}

/*
 *
 */
void c_nullAssignment(double x_data[], int x_size[2], const int idx[2])
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
  b_data[0] = false;
  b_data[1] = false;
  b_data[2] = false;
  b_data[3] = false;
  b_data[idx[0] - 1] = true;
  b_data[idx[1] - 1] = true;
  n = 0;
  i = 0;
  for (k = 0; k < 4; k++) {
    b = b_data[k];
    n += b;
    if (!b) {
      for (j = 0; j < ncolx; j++) {
        b_i = x_size[0] * j;
        x_data[i + b_i] = x_data[k + b_i];
      }
      i++;
    }
  }
  if (1 > 4 - n) {
    ncolx = 0;
  } else {
    ncolx = 4 - n;
  }
  n = x_size[1] - 1;
  for (b_i = 0; b_i <= n; b_i++) {
    for (i = 0; i < ncolx; i++) {
      x_data[i + ncolx * b_i] = x_data[i + x_size[0] * b_i];
    }
  }
  x_size[0] = ncolx;
  x_size[1] = n + 1;
}

/*
 *
 */
void nullAssignment(emxArray_real_T *x, int idx)
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

/* End of code generation (nullAssignment.c) */
