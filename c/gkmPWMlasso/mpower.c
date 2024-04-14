/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mpower.c
 *
 * Code generation for function 'mpower'
 *
 */

/* Include files */
#include "mpower.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include "cblas.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void mpower(const emxArray_real_T *a, emxArray_real_T *c)
{
  blasint idxmax_t;
  emxArray_int32_T *ipiv;
  emxArray_int32_T *p;
  emxArray_real_T *x;
  const double *a_data;
  double temp;
  double *c_data;
  double *x_data;
  int b;
  int b_n;
  int i;
  int j;
  int jj;
  int jp1j;
  int k;
  int mmj_tmp;
  int n;
  int temp_tmp;
  int u1;
  int yk;
  int *ipiv_data;
  int *p_data;
  a_data = a->data;
  if ((a->size[0] == 0) || (a->size[1] == 0)) {
    i = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity_real_T(c, i);
    c_data = c->data;
    yk = a->size[0] * a->size[1];
    for (i = 0; i < yk; i++) {
      c_data[i] = a_data[i];
    }
  } else {
    n = a->size[0];
    i = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity_real_T(c, i);
    c_data = c->data;
    yk = a->size[0] * a->size[1];
    for (i = 0; i < yk; i++) {
      c_data[i] = 0.0;
    }
    emxInit_real_T(&x, 2);
    i = x->size[0] * x->size[1];
    x->size[0] = a->size[0];
    x->size[1] = a->size[1];
    emxEnsureCapacity_real_T(x, i);
    x_data = x->data;
    yk = a->size[0] * a->size[1];
    for (i = 0; i < yk; i++) {
      x_data[i] = a_data[i];
    }
    emxInit_int32_T(&ipiv, 2);
    b_n = a->size[0];
    i = ipiv->size[0] * ipiv->size[1];
    ipiv->size[0] = 1;
    ipiv->size[1] = a->size[0];
    emxEnsureCapacity_int32_T(ipiv, i);
    ipiv_data = ipiv->data;
    ipiv_data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      ipiv_data[k - 1] = yk;
    }
    yk = a->size[0] - 1;
    u1 = a->size[0];
    if (yk <= u1) {
      u1 = yk;
    }
    for (j = 0; j < u1; j++) {
      mmj_tmp = n - j;
      b = j * (n + 1);
      jj = j * (a->size[0] + 1);
      jp1j = b + 2;
      if (mmj_tmp < 1) {
        yk = -1;
      } else {
        idxmax_t = cblas_idamax((blasint)mmj_tmp, &x_data[b], (blasint)1);
        yk = (int)idxmax_t;
      }
      if (x_data[jj + yk] != 0.0) {
        if (yk != 0) {
          b_n = j + yk;
          ipiv_data[j] = b_n + 1;
          for (k = 0; k < n; k++) {
            yk = k * n;
            temp_tmp = j + yk;
            temp = x_data[temp_tmp];
            i = b_n + yk;
            x_data[temp_tmp] = x_data[i];
            x_data[i] = temp;
          }
        }
        i = jj + mmj_tmp;
        for (yk = jp1j; yk <= i; yk++) {
          x_data[yk - 1] /= x_data[jj];
        }
      }
      if (mmj_tmp - 1 >= 1) {
        cblas_dger(CblasColMajor, (blasint)(mmj_tmp - 1),
                   (blasint)(mmj_tmp - 1), -1.0, &x_data[jj + 1], (blasint)1,
                   &x_data[b + n], (blasint)n, &x_data[(b + n) + 1],
                   (blasint)n);
      }
    }
    emxInit_int32_T(&p, 2);
    b_n = a->size[0];
    i = p->size[0] * p->size[1];
    p->size[0] = 1;
    p->size[1] = a->size[0];
    emxEnsureCapacity_int32_T(p, i);
    p_data = p->data;
    p_data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      p_data[k - 1] = yk;
    }
    i = ipiv->size[1];
    for (k = 0; k < i; k++) {
      b_n = ipiv_data[k];
      if (b_n > k + 1) {
        yk = p_data[b_n - 1];
        p_data[b_n - 1] = p_data[k];
        p_data[k] = yk;
      }
    }
    emxFree_int32_T(&ipiv);
    for (k = 0; k < n; k++) {
      i = p_data[k];
      c_data[k + c->size[0] * (i - 1)] = 1.0;
      for (j = k + 1; j <= n; j++) {
        if (c_data[(j + c->size[0] * (i - 1)) - 1] != 0.0) {
          b_n = j + 1;
          for (yk = b_n; yk <= n; yk++) {
            c_data[(yk + c->size[0] * (i - 1)) - 1] -=
                c_data[(j + c->size[0] * (i - 1)) - 1] *
                x_data[(yk + x->size[0] * (j - 1)) - 1];
          }
        }
      }
    }
    emxFree_int32_T(&p);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
                CblasNonUnit, (blasint)a->size[0], (blasint)a->size[0], 1.0,
                &x_data[0], (blasint)a->size[0], &c_data[0],
                (blasint)a->size[0]);
    emxFree_real_T(&x);
  }
}

/* End of code generation (mpower.c) */
