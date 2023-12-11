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
#include "mapTF2_ls_emxutil.h"
#include "mapTF2_ls_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_maximum(const emxArray_real_T *x, emxArray_real_T *ex)
{
  const double *x_data;
  double a;
  double b;
  double *ex_data;
  int j;
  int n;
  bool p;
  x_data = x->data;
  n = x->size[1];
  j = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = x->size[1];
  emxEnsureCapacity_real_T(ex, j);
  ex_data = ex->data;
  if (x->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      a = x_data[4 * j];
      b = x_data[4 * j + 1];
      if (rtIsNaN(b)) {
        p = false;
      } else if (rtIsNaN(a)) {
        p = true;
      } else {
        p = (a < b);
      }
      if (p) {
        a = b;
      }
      b = x_data[4 * j + 2];
      if (rtIsNaN(b)) {
        p = false;
      } else if (rtIsNaN(a)) {
        p = true;
      } else {
        p = (a < b);
      }
      if (p) {
        a = b;
      }
      b = x_data[4 * j + 3];
      if (rtIsNaN(b)) {
        p = false;
      } else if (rtIsNaN(a)) {
        p = true;
      } else {
        p = (a < b);
      }
      if (p) {
        a = b;
      }
      ex_data[j] = a;
    }
  }
}

/*
 *
 */
void c_maximum(const emxArray_real_T *x, double *ex, int *idx)
{
  const double *x_data;
  double d;
  int i;
  int k;
  int last;
  bool exitg1;
  x_data = x->data;
  last = x->size[0];
  if (x->size[0] <= 2) {
    if (x->size[0] == 1) {
      *ex = x_data[0];
      *idx = 1;
    } else if ((x_data[0] < x_data[x->size[0] - 1]) ||
               (rtIsNaN(x_data[0]) && (!rtIsNaN(x_data[x->size[0] - 1])))) {
      *ex = x_data[x->size[0] - 1];
      *idx = x->size[0];
    } else {
      *ex = x_data[0];
      *idx = 1;
    }
  } else {
    if (!rtIsNaN(x_data[0])) {
      *idx = 1;
    } else {
      *idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!rtIsNaN(x_data[k - 1])) {
          *idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (*idx == 0) {
      *ex = x_data[0];
      *idx = 1;
    } else {
      *ex = x_data[*idx - 1];
      i = *idx + 1;
      for (k = i; k <= last; k++) {
        d = x_data[k - 1];
        if (*ex < d) {
          *ex = d;
          *idx = k;
        }
      }
    }
  }
}

int d_binary_expand_op(emxArray_real_T *mat, int i1, int i,
                       const emxArray_real_T *TRANS, int i2,
                       const emxArray_real_T *C, int j,
                       const emxArray_real_T *pwm_prob)
{
  emxArray_real_T *b_mat;
  const double *C_data;
  const double *TRANS_data;
  const double *pwm_prob_data;
  double b_pwm_prob;
  double *b_mat_data;
  double *mat_data;
  int C_tmp;
  int b_i;
  int iindx;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  pwm_prob_data = pwm_prob->data;
  C_data = C->data;
  TRANS_data = TRANS->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 1);
  C_tmp = (int)(C_data[j] + 1.0);
  b_pwm_prob = pwm_prob_data[(C_tmp + pwm_prob->size[0] * (i + 1)) - 1];
  b_i = b_mat->size[0];
  if (i2 + 1 == 1) {
    b_mat->size[0] = i1 + 1;
  } else {
    b_mat->size[0] = i2 + 1;
  }
  emxEnsureCapacity_real_T(b_mat, b_i);
  b_mat_data = b_mat->data;
  stride_0_0 = (i1 + 1 != 1);
  stride_1_0 = (i2 + 1 != 1);
  if (i2 + 1 == 1) {
    loop_ub = i1 + 1;
  } else {
    loop_ub = i2 + 1;
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    b_mat_data[b_i] =
        (mat_data[b_i * stride_0_0 + mat->size[0] * i] +
         TRANS_data[b_i * stride_1_0 + TRANS->size[0] * (C_tmp - 1)]) +
        b_pwm_prob;
  }
  c_maximum(b_mat, &b_pwm_prob, &iindx);
  mat_data[((int)(C_data[j] + 1.0) + mat->size[0] * (i + 1)) - 1] = b_pwm_prob;
  emxFree_real_T(&b_mat);
  return iindx;
}

int e_binary_expand_op(emxArray_real_T *mat, int i1, int i,
                       const emxArray_real_T *TRANS, int i2,
                       const emxArray_real_T *neg, double y)
{
  emxArray_real_T *b_mat;
  const double *TRANS_data;
  const double *neg_data;
  double b_neg;
  double *b_mat_data;
  double *mat_data;
  int b_TRANS;
  int b_i;
  int iindx;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  neg_data = neg->data;
  TRANS_data = TRANS->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 1);
  b_TRANS = TRANS->size[1];
  b_neg = neg_data[i];
  b_i = b_mat->size[0];
  if (i2 + 1 == 1) {
    b_mat->size[0] = i1 + 1;
  } else {
    b_mat->size[0] = i2 + 1;
  }
  emxEnsureCapacity_real_T(b_mat, b_i);
  b_mat_data = b_mat->data;
  stride_0_0 = (i1 + 1 != 1);
  stride_1_0 = (i2 + 1 != 1);
  if (i2 + 1 == 1) {
    loop_ub = i1 + 1;
  } else {
    loop_ub = i2 + 1;
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    b_mat_data[b_i] =
        (mat_data[b_i * stride_0_0 + mat->size[0] * i] +
         TRANS_data[b_i * stride_1_0 + TRANS->size[0] * (b_TRANS - 1)]) +
        b_neg;
  }
  c_maximum(b_mat, &b_neg, &iindx);
  mat_data[((int)(y + 1.0) + mat->size[0] * (i + 1)) - 1] = b_neg;
  emxFree_real_T(&b_mat);
  return iindx;
}

/*
 *
 */
double maximum(const emxArray_real_T *x)
{
  const double *x_data;
  double d;
  double ex;
  int idx;
  int k;
  int last;
  bool exitg1;
  x_data = x->data;
  last = x->size[0];
  if (x->size[0] <= 2) {
    if (x->size[0] == 1) {
      ex = x_data[0];
    } else if ((x_data[0] < x_data[x->size[0] - 1]) ||
               (rtIsNaN(x_data[0]) && (!rtIsNaN(x_data[x->size[0] - 1])))) {
      ex = x_data[x->size[0] - 1];
    } else {
      ex = x_data[0];
    }
  } else {
    if (!rtIsNaN(x_data[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!rtIsNaN(x_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x_data[0];
    } else {
      ex = x_data[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        d = x_data[k - 1];
        if (ex < d) {
          ex = d;
        }
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
  double a;
  double b;
  double *ex_data;
  int j;
  int n;
  bool p;
  x_data = x->data;
  n = x->size[1];
  j = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = x->size[1];
  emxEnsureCapacity_real_T(ex, j);
  ex_data = ex->data;
  if (x->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      a = x_data[4 * j];
      b = x_data[4 * j + 1];
      if (rtIsNaN(b)) {
        p = false;
      } else if (rtIsNaN(a)) {
        p = true;
      } else {
        p = (a > b);
      }
      if (p) {
        a = b;
      }
      b = x_data[4 * j + 2];
      if (rtIsNaN(b)) {
        p = false;
      } else if (rtIsNaN(a)) {
        p = true;
      } else {
        p = (a > b);
      }
      if (p) {
        a = b;
      }
      b = x_data[4 * j + 3];
      if (rtIsNaN(b)) {
        p = false;
      } else if (rtIsNaN(a)) {
        p = true;
      } else {
        p = (a > b);
      }
      if (p) {
        a = b;
      }
      ex_data[j] = a;
    }
  }
}

/* End of code generation (minOrMax.c) */
