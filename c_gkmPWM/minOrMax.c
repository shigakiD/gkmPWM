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
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "relop.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_maximum(const double x[9], double *ex, int *idx)
{
  double d;
  int k;
  *ex = x[0];
  *idx = 1;
  for (k = 0; k < 8; k++) {
    d = x[k + 1];
    if (*ex < d) {
      *ex = d;
      *idx = k + 2;
    }
  }
}

/*
 *
 */
double b_minimum(const double x[2])
{
  double ex;
  if (x[0] > x[1]) {
    ex = x[1];
  } else {
    ex = x[0];
  }
  return ex;
}

/*
 *
 */
void c_maximum(const double x[2], double *ex, int *idx)
{
  if (x[0] < x[1]) {
    *ex = x[1];
    *idx = 2;
  } else {
    *ex = x[0];
    *idx = 1;
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
double e_maximum(const double x[4])
{
  double ex;
  ex = x[0];
  if (x[0] < x[1]) {
    ex = x[1];
  }
  if (ex < x[2]) {
    ex = x[2];
  }
  if (ex < x[3]) {
    ex = x[3];
  }
  return ex;
}

/*
 *
 */
void f_maximum(const emxArray_real_T *x, emxArray_real_T *ex)
{
  const double *x_data;
  double d;
  double d1;
  double *ex_data;
  int j;
  int n;
  x_data = x->data;
  n = x->size[1];
  j = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = x->size[1];
  emxEnsureCapacity_real_T(ex, j);
  ex_data = ex->data;
  if (x->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      d = x_data[4 * j];
      d1 = x_data[4 * j + 1];
      if (d < d1) {
        d = d1;
      }
      d1 = x_data[4 * j + 2];
      if (d < d1) {
        d = d1;
      }
      d1 = x_data[4 * j + 3];
      if (d < d1) {
        d = d1;
      }
      ex_data[j] = d;
    }
  }
}

/*
 *
 */
void g_maximum(const emxArray_real_T *x, emxArray_real_T *ex)
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
double h_maximum(const emxArray_real_T *x)
{
  const double *x_data;
  double d;
  double ex;
  int k;
  int last;
  x_data = x->data;
  last = x->size[1];
  if (x->size[1] <= 2) {
    if (x->size[1] == 1) {
      ex = x_data[0];
    } else if (x_data[0] < x_data[x->size[1] - 1]) {
      ex = x_data[x->size[1] - 1];
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
double i_maximum(const double x[9])
{
  double d;
  double ex;
  int k;
  ex = x[0];
  for (k = 0; k < 8; k++) {
    d = x[k + 1];
    if (ex < d) {
      ex = d;
    }
  }
  return ex;
}

/*
 *
 */
void j_maximum(const creal_T x_data[], creal_T *ex, int *idx)
{
  double x;
  double y;
  *idx = 1;
  *ex = x_data[0];
  absRelopProxies(x_data[0], x_data[1], &x, &y);
  if (x < y) {
    *ex = x_data[1];
    *idx = 2;
  }
  absRelopProxies(*ex, x_data[2], &x, &y);
  if (x < y) {
    *ex = x_data[2];
    *idx = 3;
  }
  absRelopProxies(*ex, x_data[3], &x, &y);
  if (x < y) {
    *ex = x_data[3];
    *idx = 4;
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
double minimum(const emxArray_real_T *x)
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
    } else if (x_data[0] > x_data[x->size[0] - 1]) {
      ex = x_data[x->size[0] - 1];
    } else {
      ex = x_data[0];
    }
  } else {
    ex = x_data[0];
    for (k = 2; k <= last; k++) {
      d = x_data[k - 1];
      if (ex > d) {
        ex = d;
      }
    }
  }
  return ex;
}

/* End of code generation (minOrMax.c) */
