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
#include "gkmPWM_rtwutil.h"
#include "gkmPWM_types.h"
#include <math.h>
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
creal_T c_minimum(const creal_T x[4])
{
  creal_T dc;
  creal_T ex;
  double absx;
  double b_x;
  double bi;
  int b_exponent;
  int c_exponent;
  int exponent;
  int k;
  bool SCALEA;
  bool SCALEB;
  ex = x[0];
  for (k = 0; k < 3; k++) {
    if ((fabs(ex.re) > 8.9884656743115785E+307) ||
        (fabs(ex.im) > 8.9884656743115785E+307)) {
      SCALEA = true;
    } else {
      SCALEA = false;
    }
    dc = x[k + 1];
    if ((fabs(dc.re) > 8.9884656743115785E+307) ||
        (fabs(x[k + 1].im) > 8.9884656743115785E+307)) {
      SCALEB = true;
    } else {
      SCALEB = false;
    }
    if (SCALEA || SCALEB) {
      b_x = rt_hypotd(ex.re / 2.0, ex.im / 2.0);
      bi = rt_hypotd(x[k + 1].re / 2.0, x[k + 1].im / 2.0);
    } else {
      b_x = rt_hypotd(ex.re, ex.im);
      bi = rt_hypotd(x[k + 1].re, x[k + 1].im);
    }
    absx = bi / 2.0;
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
    if (fabs(bi - b_x) < absx) {
      b_x = atan2(ex.im, ex.re);
      bi = atan2(x[k + 1].im, x[k + 1].re);
      absx = fabs(bi / 2.0);
      if (absx <= 2.2250738585072014E-308) {
        absx = 4.94065645841247E-324;
      } else {
        frexp(absx, &b_exponent);
        absx = ldexp(1.0, b_exponent - 53);
      }
      if (fabs(bi - b_x) < absx) {
        absx = x[k + 1].re;
        bi = x[k + 1].im;
        if (ex.re != absx) {
          if (b_x >= 0.0) {
            b_x = absx;
            bi = ex.re;
          } else {
            b_x = ex.re;
            bi = absx;
          }
        } else if (ex.re < 0.0) {
          b_x = bi;
          bi = ex.im;
        } else {
          b_x = ex.im;
        }
        absx = fabs(bi / 2.0);
        if (absx <= 2.2250738585072014E-308) {
          absx = 4.94065645841247E-324;
        } else {
          frexp(absx, &c_exponent);
          absx = ldexp(1.0, c_exponent - 53);
        }
        if (fabs(bi - b_x) < absx) {
          b_x = 0.0;
          bi = 0.0;
        }
      }
    }
    if (b_x > bi) {
      ex = dc;
    }
  }
  return ex;
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
  double d;
  double ex;
  int k;
  ex = x[0];
  for (k = 0; k < 3; k++) {
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
void f_maximum(const emxArray_real_T *x, emxArray_real_T *ex)
{
  const double *x_data;
  double d;
  double d1;
  double *ex_data;
  int i;
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
      for (i = 0; i < 3; i++) {
        d1 = x_data[(i + 4 * j) + 1];
        if (d < d1) {
          d = d1;
        }
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

void q_binary_expand_op(emxArray_real_T *Rd, const emxArray_real_T *CT,
                        const emxArray_real_T *ct)
{
  emxArray_real_T *b_CT;
  const double *CT_data;
  const double *ct_data;
  double *b_CT_data;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  ct_data = ct->data;
  CT_data = CT->data;
  emxInit_real_T(&b_CT, 2);
  i = b_CT->size[0] * b_CT->size[1];
  if (ct->size[0] == 1) {
    b_CT->size[0] = CT->size[0];
  } else {
    b_CT->size[0] = ct->size[0];
  }
  if (ct->size[1] == 1) {
    b_CT->size[1] = CT->size[1];
  } else {
    b_CT->size[1] = ct->size[1];
  }
  emxEnsureCapacity_real_T(b_CT, i);
  b_CT_data = b_CT->data;
  stride_0_0 = (CT->size[0] != 1);
  stride_0_1 = (CT->size[1] != 1);
  stride_1_0 = (ct->size[0] != 1);
  stride_1_1 = (ct->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (ct->size[1] == 1) {
    loop_ub = CT->size[1];
  } else {
    loop_ub = ct->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    if (ct->size[0] == 1) {
      b_loop_ub = CT->size[0];
    } else {
      b_loop_ub = ct->size[0];
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_CT_data[i1 + b_CT->size[0] * i] =
          CT_data[i1 * stride_0_0 + CT->size[0] * aux_0_1] -
          ct_data[i1 * stride_1_0 + ct->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  g_maximum(b_CT, Rd);
  emxFree_real_T(&b_CT);
}

/* End of code generation (minOrMax.c) */
