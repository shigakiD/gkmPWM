/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minOrMax.c
 *
 * Code generation for function 'minOrMax'
 *
 */

/* Include files */
#include "minOrMax.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_maximum(const emxArray_real_T *x, double *ex, int *idx)
{
  const double *x_data;
  double d;
  int k;
  int last;
  x_data = x->data;
  last = x->size[1];
  if (x->size[1] <= 2) {
    if (x->size[1] == 1) {
      *ex = x_data[0];
      *idx = 1;
    } else if (x_data[0] < x_data[x->size[1] - 1]) {
      *ex = x_data[x->size[1] - 1];
      *idx = x->size[1];
    } else {
      *ex = x_data[0];
      *idx = 1;
    }
  } else {
    *ex = x_data[0];
    *idx = 1;
    for (k = 2; k <= last; k++) {
      d = x_data[k - 1];
      if (*ex < d) {
        *ex = d;
        *idx = k;
      }
    }
  }
}

/*
 *
 */
void c_maximum(const emxArray_real_T *x, emxArray_real_T *ex)
{
  const double *x_data;
  double d;
  double *ex_data;
  int i;
  int j;
  int m;
  int n;
  x_data = x->data;
  m = x->size[0];
  n = x->size[1];
  j = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = x->size[1];
  emxEnsureCapacity_real_T(ex, j);
  ex_data = ex->data;
  if (x->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      ex_data[j] = x_data[x->size[0] * j];
      for (i = 2; i <= m; i++) {
        d = x_data[(i + x->size[0] * j) - 1];
        if (ex_data[j] < d) {
          ex_data[j] = d;
        }
      }
    }
  }
}

/*
 *
 */
void d_maximum(const emxArray_real_T *x, double *ex, int *idx)
{
  const double *x_data;
  double d;
  int k;
  int last;
  x_data = x->data;
  last = x->size[0];
  if (x->size[0] <= 2) {
    if (x->size[0] == 1) {
      *ex = x_data[0];
      *idx = 1;
    } else if (x_data[0] < x_data[x->size[0] - 1]) {
      *ex = x_data[x->size[0] - 1];
      *idx = x->size[0];
    } else {
      *ex = x_data[0];
      *idx = 1;
    }
  } else {
    *ex = x_data[0];
    *idx = 1;
    for (k = 2; k <= last; k++) {
      d = x_data[k - 1];
      if (*ex < d) {
        *ex = d;
        *idx = k;
      }
    }
  }
}

/*
 *
 */
double maximum(const emxArray_real_T *x)
{
  const double *x_data;
  double d;
  double ex;
  int k;
  int last;
  x_data = x->data;
  last = x->size[0];
  if (x->size[0] <= 2) {
    if (x->size[0] == 1) {
      ex = x_data[0];
    } else if (x_data[0] < x_data[x->size[0] - 1]) {
      ex = x_data[x->size[0] - 1];
    } else {
      ex = x_data[0];
    }
  } else {
    ex = x_data[0];
    for (k = 2; k <= last; k++) {
      d = x_data[k - 1];
      if (ex < d) {
        ex = d;
      }
    }
  }
  return ex;
}

/*
 *
 */
void minimum(const emxArray_real_T *x, emxArray_real_T *ex)
{
  const double *x_data;
  double d;
  double *ex_data;
  int i;
  int j;
  int m;
  int n;
  x_data = x->data;
  m = x->size[0];
  n = x->size[1];
  j = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = x->size[1];
  emxEnsureCapacity_real_T(ex, j);
  ex_data = ex->data;
  if (x->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      ex_data[j] = x_data[x->size[0] * j];
      for (i = 2; i <= m; i++) {
        d = x_data[(i + x->size[0] * j) - 1];
        if (ex_data[j] > d) {
          ex_data[j] = d;
        }
      }
    }
  }
}

/* End of code generation (minOrMax.c) */
