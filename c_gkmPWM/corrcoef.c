/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * corrcoef.c
 *
 * Code generation for function 'corrcoef'
 *
 */

/* Include files */
#include "corrcoef.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_corrcoef(const emxArray_real_T *x, emxArray_real_T *r)
{
  emxArray_real_T *b_x;
  emxArray_real_T *d;
  const double *x_data;
  double b_d;
  double s;
  double *d_data;
  double *r_data;
  int b_i;
  int b_n;
  int fm;
  int i;
  int j;
  int k;
  int n;
  x_data = x->data;
  emxInit_real_T(&b_x, 2);
  i = b_x->size[0] * b_x->size[1];
  b_x->size[0] = x->size[0];
  b_x->size[1] = x->size[1];
  emxEnsureCapacity_real_T(b_x, i);
  d_data = b_x->data;
  fm = x->size[0] * x->size[1];
  for (i = 0; i < fm; i++) {
    d_data[i] = x_data[i];
  }
  n = x->size[0] - 1;
  b_n = x->size[1];
  i = r->size[0] * r->size[1];
  r->size[0] = x->size[1];
  r->size[1] = x->size[1];
  emxEnsureCapacity_real_T(r, i);
  r_data = r->data;
  fm = x->size[1] * x->size[1];
  for (i = 0; i < fm; i++) {
    r_data[i] = 0.0;
  }
  if (x->size[0] < 2) {
    i = r->size[0] * r->size[1];
    r->size[0] = x->size[1];
    r->size[1] = x->size[1];
    emxEnsureCapacity_real_T(r, i);
    r_data = r->data;
    fm = x->size[1] * x->size[1];
    for (i = 0; i < fm; i++) {
      r_data[i] = 0.0;
    }
  } else {
    for (j = 0; j < b_n; j++) {
      s = 0.0;
      for (b_i = 0; b_i <= n; b_i++) {
        s += d_data[b_i + b_x->size[0] * j];
      }
      s /= (double)x->size[0];
      for (b_i = 0; b_i <= n; b_i++) {
        d_data[b_i + b_x->size[0] * j] -= s;
      }
    }
    fm = x->size[0] - 1;
    for (j = 0; j < b_n; j++) {
      s = 0.0;
      for (k = 0; k <= n; k++) {
        b_d = d_data[k + b_x->size[0] * j];
        s += b_d * b_d;
      }
      r_data[j + r->size[0] * j] = s / (double)fm;
      i = j + 2;
      for (b_i = i; b_i <= b_n; b_i++) {
        s = 0.0;
        for (k = 0; k <= n; k++) {
          s += d_data[k + b_x->size[0] * (b_i - 1)] *
               d_data[k + b_x->size[0] * j];
        }
        r_data[(b_i + r->size[0] * j) - 1] = s / (double)fm;
      }
    }
  }
  emxFree_real_T(&b_x);
  emxInit_real_T(&d, 1);
  fm = r->size[0];
  i = d->size[0];
  d->size[0] = r->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  for (k = 0; k < fm; k++) {
    d_data[k] = sqrt(r_data[k + r->size[0] * k]);
  }
  for (j = 0; j < fm; j++) {
    i = j + 2;
    for (b_i = i; b_i <= fm; b_i++) {
      r_data[(b_i + r->size[0] * j) - 1] =
          r_data[(b_i + r->size[0] * j) - 1] / d_data[b_i - 1] / d_data[j];
    }
    for (b_i = i; b_i <= fm; b_i++) {
      b_d = r_data[(b_i + r->size[0] * j) - 1];
      s = fabs(b_d);
      if (s > 1.0) {
        r_data[(b_i + r->size[0] * j) - 1] = b_d / s;
      }
      r_data[j + r->size[0] * (b_i - 1)] = r_data[(b_i + r->size[0] * j) - 1];
    }
    if (r_data[j + r->size[0] * j] > 0.0) {
      s = r_data[j + r->size[0] * j];
      if (s < 0.0) {
        s = -1.0;
      } else if (s > 0.0) {
        s = 1.0;
      }
      r_data[j + r->size[0] * j] = s;
    } else {
      r_data[j + r->size[0] * j] = 0.0;
    }
  }
  emxFree_real_T(&d);
}

/*
 *
 */
void corrcoef(const emxArray_real_T *x, const emxArray_real_T *varargin_1,
              double r[4])
{
  emxArray_real_T *result;
  double d[2];
  const double *varargin_1_data;
  const double *x_data;
  double b_d;
  double s;
  double *result_data;
  int b_i;
  int fm;
  int i;
  int j;
  int k;
  int n;
  int r_tmp;
  varargin_1_data = varargin_1->data;
  x_data = x->data;
  emxInit_real_T(&result, 2);
  i = result->size[0] * result->size[1];
  result->size[0] = x->size[0];
  result->size[1] = 2;
  emxEnsureCapacity_real_T(result, i);
  result_data = result->data;
  fm = x->size[0];
  for (i = 0; i < fm; i++) {
    result_data[i] = x_data[i];
  }
  fm = varargin_1->size[0];
  for (i = 0; i < fm; i++) {
    result_data[i + result->size[0]] = varargin_1_data[i];
  }
  fm = result->size[0];
  n = result->size[0] - 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  if (result->size[0] >= 2) {
    for (j = 0; j < 2; j++) {
      s = 0.0;
      for (b_i = 0; b_i <= n; b_i++) {
        s += result_data[b_i + result->size[0] * j];
      }
      s /= (double)fm;
      for (b_i = 0; b_i <= n; b_i++) {
        result_data[b_i + result->size[0] * j] -= s;
      }
    }
    fm--;
    for (j = 0; j < 2; j++) {
      s = 0.0;
      for (k = 0; k <= n; k++) {
        b_d = result_data[k + result->size[0] * j];
        s += b_d * b_d;
      }
      r_tmp = j << 1;
      r[j + r_tmp] = s / (double)fm;
      i = j + 2;
      for (b_i = i; b_i < 3; b_i++) {
        s = 0.0;
        for (k = 0; k <= n; k++) {
          s += result_data[k + result->size[0]] *
               result_data[k + result->size[0] * j];
        }
        r[r_tmp + 1] = s / (double)fm;
      }
    }
  }
  emxFree_real_T(&result);
  d[0] = sqrt(r[0]);
  d[1] = sqrt(r[3]);
  for (j = 0; j < 2; j++) {
    i = j + 2;
    for (b_i = i; b_i < 3; b_i++) {
      r_tmp = (j << 1) + 1;
      r[r_tmp] = r[r_tmp] / d[1] / d[j];
    }
    for (b_i = i; b_i < 3; b_i++) {
      fm = (j << 1) + 1;
      b_d = r[fm];
      s = fabs(b_d);
      if (s > 1.0) {
        r[fm] = b_d / s;
      }
      r[j + 2] = r[fm];
    }
    i = j + (j << 1);
    r[i] = (r[i] > 0.0);
  }
}

/* End of code generation (corrcoef.c) */
