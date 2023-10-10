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
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void corrcoef(const emxArray_real_T *x, const emxArray_real_T *varargin_1,
              double b_r[4])
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
  b_r[0] = 0.0;
  b_r[1] = 0.0;
  b_r[2] = 0.0;
  b_r[3] = 0.0;
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
      b_r[j + r_tmp] = s / (double)fm;
      i = j + 2;
      for (b_i = i; b_i < 3; b_i++) {
        s = 0.0;
        for (k = 0; k <= n; k++) {
          s += result_data[k + result->size[0]] *
               result_data[k + result->size[0] * j];
        }
        b_r[r_tmp + 1] = s / (double)fm;
      }
    }
  }
  emxFree_real_T(&result);
  d[0] = sqrt(b_r[0]);
  d[1] = sqrt(b_r[3]);
  for (j = 0; j < 2; j++) {
    i = j + 2;
    for (b_i = i; b_i < 3; b_i++) {
      r_tmp = (j << 1) + 1;
      b_r[r_tmp] = b_r[r_tmp] / d[1] / d[j];
    }
    for (b_i = i; b_i < 3; b_i++) {
      fm = (j << 1) + 1;
      b_d = b_r[fm];
      s = fabs(b_d);
      if (s > 1.0) {
        b_r[fm] = b_d / s;
      }
      b_r[j + 2] = b_r[fm];
    }
    i = j + (j << 1);
    b_r[i] = (b_r[i] > 0.0);
  }
}

/* End of code generation (corrcoef.c) */
