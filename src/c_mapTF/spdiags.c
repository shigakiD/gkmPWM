/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spdiags.c
 *
 * Code generation for function 'spdiags'
 *
 */

/* Include files */
#include "spdiags.h"
#include "colon.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include "unique.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void diagk(const emxArray_real_T *X, double k, emxArray_real_T *D);

/* Function Definitions */
/*
 *
 */
static void diagk(const emxArray_real_T *X, double k, emxArray_real_T *D)
{
  const double *X_data;
  double *D_data;
  int absk;
  int i;
  int u0;
  int u1;
  bool negk;
  X_data = X->data;
  if ((X->size[0] != 1) && (X->size[1] != 1)) {
    if (k < 0.0) {
      absk = (int)-k;
      negk = true;
    } else {
      absk = (int)k;
      negk = false;
    }
    if ((X->size[0] == 1) && (X->size[1] == 1) && (absk == 0)) {
      u0 = D->size[0];
      D->size[0] = 1;
      emxEnsureCapacity_real_T(D, u0);
      D_data = D->data;
      D_data[0] = X_data[0];
    } else if ((negk && (absk > X->size[0])) ||
               ((!negk) && (absk > X->size[1]))) {
      D->size[0] = 0;
    } else {
      if (negk) {
        if (absk < X->size[0]) {
          u0 = X->size[0] - absk;
          u1 = X->size[1];
          if (u0 <= u1) {
            u1 = u0;
          }
          i = absk;
          absk = 0;
        } else {
          u1 = 0;
          i = 0;
          absk = 0;
        }
      } else if (absk < X->size[1]) {
        u0 = X->size[0];
        u1 = X->size[1] - absk;
        if (u0 <= u1) {
          u1 = u0;
        }
        i = 0;
      } else {
        u1 = 0;
        i = 0;
        absk = 0;
      }
      u0 = D->size[0];
      D->size[0] = u1;
      emxEnsureCapacity_real_T(D, u0);
      D_data = D->data;
      u0 = u1 - 1;
      for (u1 = 0; u1 <= u0; u1++) {
        D_data[u1] = X_data[(i + u1) + X->size[0] * (absk + u1)];
      }
    }
  } else if ((X->size[0] != 0) && (X->size[1] != 0) && (0.0 <= k) &&
             (k + 1.0 <= X->size[1])) {
    u0 = D->size[0];
    D->size[0] = 1;
    emxEnsureCapacity_real_T(D, u0);
    D_data = D->data;
    D_data[0] = X_data[(int)(k + 1.0) - 1];
  } else if ((X->size[0] != 0) && (X->size[1] != 0) && (k < 0.0) &&
             (1.0 - k <= X->size[0])) {
    u0 = D->size[0];
    D->size[0] = 1;
    emxEnsureCapacity_real_T(D, u0);
    D_data = D->data;
    D_data[0] = X_data[(int)(1.0 - k) - 1];
  } else {
    D->size[0] = 0;
  }
}

/*
 *
 */
void spdiags(const emxArray_real_T *arg1, emxArray_real_T *res1)
{
  emxArray_int32_T *i_tmp;
  emxArray_int32_T *j;
  emxArray_real_T *b_j;
  emxArray_real_T *d;
  emxArray_real_T *i;
  const double *arg1_data;
  double a;
  double b;
  double *b_j_data;
  double *d_data;
  double *i_data;
  double *res1_data;
  int b_i;
  int idx;
  int ii;
  int jj;
  int m;
  int n;
  int nx;
  int *i_tmp_data;
  int *j_data;
  bool exitg1;
  bool guard1 = false;
  arg1_data = arg1->data;
  emxInit_real_T(&d, 1);
  emxInit_real_T(&i, 2);
  i_data = i->data;
  emxInit_int32_T(&j, 1);
  j_data = j->data;
  emxInit_int32_T(&i_tmp, 1);
  i_tmp_data = i_tmp->data;
  emxInit_real_T(&b_j, 1);
  if ((arg1->size[0] == 0) || (arg1->size[1] == 0)) {
    res1->size[0] = 0;
    res1->size[1] = 0;
  } else {
    n = arg1->size[1];
    m = arg1->size[0];
    nx = arg1->size[0] * arg1->size[1];
    if (nx == 0) {
      i_tmp->size[0] = 0;
      j->size[0] = 0;
    } else {
      idx = 0;
      b_i = i_tmp->size[0];
      i_tmp->size[0] = nx;
      emxEnsureCapacity_int32_T(i_tmp, b_i);
      i_tmp_data = i_tmp->data;
      b_i = j->size[0];
      j->size[0] = nx;
      emxEnsureCapacity_int32_T(j, b_i);
      j_data = j->data;
      ii = 1;
      jj = 1;
      exitg1 = false;
      while ((!exitg1) && (jj <= arg1->size[1])) {
        guard1 = false;
        if (arg1_data[(ii + arg1->size[0] * (jj - 1)) - 1] != 0.0) {
          idx++;
          i_tmp_data[idx - 1] = ii;
          j_data[idx - 1] = jj;
          if (idx >= nx) {
            exitg1 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
        if (guard1) {
          ii++;
          if (ii > arg1->size[0]) {
            ii = 1;
            jj++;
          }
        }
      }
      if (nx == 1) {
        if (idx == 0) {
          i_tmp->size[0] = 0;
          j->size[0] = 0;
        }
      } else {
        b_i = i_tmp->size[0];
        if (1 > idx) {
          i_tmp->size[0] = 0;
        } else {
          i_tmp->size[0] = idx;
        }
        emxEnsureCapacity_int32_T(i_tmp, b_i);
        i_tmp_data = i_tmp->data;
        b_i = j->size[0];
        if (1 > idx) {
          j->size[0] = 0;
        } else {
          j->size[0] = idx;
        }
        emxEnsureCapacity_int32_T(j, b_i);
        j_data = j->data;
      }
    }
    b_i = b_j->size[0];
    b_j->size[0] = j->size[0];
    emxEnsureCapacity_real_T(b_j, b_i);
    b_j_data = b_j->data;
    ii = j->size[0];
    for (b_i = 0; b_i < ii; b_i++) {
      b_j_data[b_i] = (double)j_data[b_i] - (double)i_tmp_data[b_i];
    }
    unique_vector(b_j, d);
    d_data = d->data;
    if (d->size[0] == 0) {
      res1->size[0] = arg1->size[0];
      res1->size[1] = 0;
    } else {
      b_i = res1->size[0] * res1->size[1];
      res1->size[0] = (int)fmin(arg1->size[0], arg1->size[1]);
      res1->size[1] = d->size[0];
      emxEnsureCapacity_real_T(res1, b_i);
      res1_data = res1->data;
      ii = (int)fmin(arg1->size[0], arg1->size[1]) * d->size[0];
      for (b_i = 0; b_i < ii; b_i++) {
        res1_data[b_i] = 0.0;
      }
      if (arg1->size[0] >= arg1->size[1]) {
        b_i = d->size[0];
        for (jj = 0; jj < b_i; jj++) {
          a = fmax(1.0, d_data[jj] + 1.0);
          b = fmin(n, (double)m + d_data[jj]);
          if (b < a) {
            i->size[0] = 1;
            i->size[1] = 0;
          } else if (floor(a) == a) {
            nx = i->size[0] * i->size[1];
            i->size[0] = 1;
            ii = (int)floor(b - a);
            i->size[1] = ii + 1;
            emxEnsureCapacity_real_T(i, nx);
            i_data = i->data;
            for (nx = 0; nx <= ii; nx++) {
              i_data[nx] = a + (double)nx;
            }
          } else {
            eml_float_colon(a, b, i);
            i_data = i->data;
          }
          diagk(arg1, d_data[jj], b_j);
          b_j_data = b_j->data;
          ii = b_j->size[0];
          for (nx = 0; nx < ii; nx++) {
            res1_data[((int)i_data[nx] + res1->size[0] * jj) - 1] =
                b_j_data[nx];
          }
        }
      } else {
        b_i = d->size[0];
        for (jj = 0; jj < b_i; jj++) {
          a = fmax(1.0, 1.0 - d_data[jj]);
          b = fmin(m, (double)n - d_data[jj]);
          if (b < a) {
            i->size[0] = 1;
            i->size[1] = 0;
          } else if (floor(a) == a) {
            nx = i->size[0] * i->size[1];
            i->size[0] = 1;
            ii = (int)floor(b - a);
            i->size[1] = ii + 1;
            emxEnsureCapacity_real_T(i, nx);
            i_data = i->data;
            for (nx = 0; nx <= ii; nx++) {
              i_data[nx] = a + (double)nx;
            }
          } else {
            eml_float_colon(a, b, i);
            i_data = i->data;
          }
          diagk(arg1, d_data[jj], b_j);
          b_j_data = b_j->data;
          ii = b_j->size[0];
          for (nx = 0; nx < ii; nx++) {
            res1_data[((int)i_data[nx] + res1->size[0] * jj) - 1] =
                b_j_data[nx];
          }
        }
      }
    }
  }
  emxFree_real_T(&b_j);
  emxFree_int32_T(&i_tmp);
  emxFree_int32_T(&j);
  emxFree_real_T(&i);
  emxFree_real_T(&d);
}

/* End of code generation (spdiags.c) */
