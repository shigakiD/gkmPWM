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
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include "cblas.h"
#include "lapacke.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void mpower(const emxArray_real_T *a, emxArray_real_T *c)
{
  lapack_int info_t;
  lapack_int *ipiv_t_data;
  emxArray_int32_T *ipiv;
  emxArray_int32_T *p;
  emxArray_lapack_int *ipiv_t;
  emxArray_real_T *x;
  const double *a_data;
  double *c_data;
  double *x_data;
  int b_i;
  int b_n;
  int i;
  int k;
  int n;
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
    emxInit_lapack_int(&ipiv_t);
    i = ipiv_t->size[0];
    ipiv_t->size[0] = a->size[0];
    emxEnsureCapacity_lapack_int(ipiv_t, i);
    ipiv_t_data = ipiv_t->data;
    info_t = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, (lapack_int)a->size[0],
                                 (lapack_int)a->size[0], &x_data[0],
                                 (lapack_int)a->size[0], &ipiv_t_data[0]);
    i = ipiv->size[0] * ipiv->size[1];
    ipiv->size[0] = 1;
    ipiv->size[1] = ipiv_t->size[0];
    emxEnsureCapacity_int32_T(ipiv, i);
    ipiv_data = ipiv->data;
    if ((int)info_t < 0) {
      yk = x->size[0];
      b_n = x->size[1];
      i = x->size[0] * x->size[1];
      x->size[0] = yk;
      x->size[1] = b_n;
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      yk *= b_n;
      for (i = 0; i < yk; i++) {
        x_data[i] = 0.0;
      }
      i = ipiv_t->size[0] - 1;
      for (k = 0; k <= i; k++) {
        ipiv_data[k] = k + 1;
      }
    } else {
      i = ipiv_t->size[0] - 1;
      for (k = 0; k <= i; k++) {
        ipiv_data[k] = (int)ipiv_t_data[k];
      }
    }
    emxFree_lapack_int(&ipiv_t);
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
      for (yk = k + 1; yk <= n; yk++) {
        if (c_data[(yk + c->size[0] * (i - 1)) - 1] != 0.0) {
          b_n = yk + 1;
          for (b_i = b_n; b_i <= n; b_i++) {
            c_data[(b_i + c->size[0] * (i - 1)) - 1] -=
                c_data[(yk + c->size[0] * (i - 1)) - 1] *
                x_data[(b_i + x->size[0] * (yk - 1)) - 1];
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
