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
  emxArray_int32_T *idx;
  emxArray_int32_T *iwork;
  emxArray_int32_T *j;
  emxArray_real_T *a;
  emxArray_real_T *d;
  emxArray_real_T *i;
  const double *arg1_data;
  double absx;
  double x;
  double *a_data;
  double *d_data;
  double *i_data;
  double *res1_data;
  int b_i;
  int b_idx;
  int b_n;
  int exitg2;
  int exponent;
  int i1;
  int ii;
  int jj;
  int k;
  int kEnd;
  int m;
  int n;
  int na;
  int nx;
  int p;
  int q;
  int qEnd;
  int *idx_data;
  int *iwork_data;
  bool exitg1;
  bool guard1 = false;
  arg1_data = arg1->data;
  emxInit_real_T(&d, 1);
  emxInit_real_T(&i, 2);
  i_data = i->data;
  emxInit_int32_T(&j, 1);
  iwork_data = j->data;
  emxInit_real_T(&a, 1);
  emxInit_int32_T(&idx, 1);
  idx_data = idx->data;
  emxInit_int32_T(&iwork, 1);
  if ((arg1->size[0] == 0) || (arg1->size[1] == 0)) {
    res1->size[0] = 0;
    res1->size[1] = 0;
  } else {
    n = arg1->size[1];
    m = arg1->size[0];
    nx = arg1->size[0] * arg1->size[1];
    if (nx == 0) {
      idx->size[0] = 0;
      j->size[0] = 0;
    } else {
      b_idx = 0;
      b_i = idx->size[0];
      idx->size[0] = nx;
      emxEnsureCapacity_int32_T(idx, b_i);
      idx_data = idx->data;
      b_i = j->size[0];
      j->size[0] = nx;
      emxEnsureCapacity_int32_T(j, b_i);
      iwork_data = j->data;
      ii = 1;
      jj = 1;
      exitg1 = false;
      while ((!exitg1) && (jj <= arg1->size[1])) {
        guard1 = false;
        if (arg1_data[(ii + arg1->size[0] * (jj - 1)) - 1] != 0.0) {
          b_idx++;
          idx_data[b_idx - 1] = ii;
          iwork_data[b_idx - 1] = jj;
          if (b_idx >= nx) {
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
        if (b_idx == 0) {
          idx->size[0] = 0;
          j->size[0] = 0;
        }
      } else {
        b_i = idx->size[0];
        if (1 > b_idx) {
          idx->size[0] = 0;
        } else {
          idx->size[0] = b_idx;
        }
        emxEnsureCapacity_int32_T(idx, b_i);
        idx_data = idx->data;
        b_i = j->size[0];
        if (1 > b_idx) {
          j->size[0] = 0;
        } else {
          j->size[0] = b_idx;
        }
        emxEnsureCapacity_int32_T(j, b_i);
        iwork_data = j->data;
      }
    }
    b_i = a->size[0];
    a->size[0] = j->size[0];
    emxEnsureCapacity_real_T(a, b_i);
    a_data = a->data;
    ii = j->size[0];
    for (b_i = 0; b_i < ii; b_i++) {
      a_data[b_i] = (double)iwork_data[b_i] - (double)idx_data[b_i];
    }
    na = a->size[0];
    b_n = a->size[0] + 1;
    b_i = idx->size[0];
    idx->size[0] = a->size[0];
    emxEnsureCapacity_int32_T(idx, b_i);
    idx_data = idx->data;
    ii = a->size[0];
    for (b_i = 0; b_i < ii; b_i++) {
      idx_data[b_i] = 0;
    }
    if (a->size[0] != 0) {
      b_i = iwork->size[0];
      iwork->size[0] = a->size[0];
      emxEnsureCapacity_int32_T(iwork, b_i);
      iwork_data = iwork->data;
      b_i = a->size[0] - 1;
      for (k = 1; k <= b_i; k += 2) {
        if (a_data[k - 1] <= a_data[k]) {
          idx_data[k - 1] = k;
          idx_data[k] = k + 1;
        } else {
          idx_data[k - 1] = k + 1;
          idx_data[k] = k;
        }
      }
      if ((a->size[0] & 1) != 0) {
        idx_data[a->size[0] - 1] = a->size[0];
      }
      ii = 2;
      while (ii < b_n - 1) {
        jj = ii << 1;
        nx = 1;
        for (b_idx = ii + 1; b_idx < b_n; b_idx = qEnd + ii) {
          p = nx;
          q = b_idx;
          qEnd = nx + jj;
          if (qEnd > b_n) {
            qEnd = b_n;
          }
          k = 0;
          kEnd = qEnd - nx;
          while (k + 1 <= kEnd) {
            b_i = idx_data[q - 1];
            i1 = idx_data[p - 1];
            if (a_data[i1 - 1] <= a_data[b_i - 1]) {
              iwork_data[k] = i1;
              p++;
              if (p == b_idx) {
                while (q < qEnd) {
                  k++;
                  iwork_data[k] = idx_data[q - 1];
                  q++;
                }
              }
            } else {
              iwork_data[k] = b_i;
              q++;
              if (q == qEnd) {
                while (p < b_idx) {
                  k++;
                  iwork_data[k] = idx_data[p - 1];
                  p++;
                }
              }
            }
            k++;
          }
          for (k = 0; k < kEnd; k++) {
            idx_data[(nx + k) - 1] = iwork_data[k];
          }
          nx = qEnd;
        }
        ii = jj;
      }
    }
    b_i = d->size[0];
    d->size[0] = a->size[0];
    emxEnsureCapacity_real_T(d, b_i);
    d_data = d->data;
    for (k = 0; k < na; k++) {
      d_data[k] = a_data[idx_data[k] - 1];
    }
    k = a->size[0];
    ii = -1;
    jj = 0;
    while (jj + 1 <= k) {
      x = d_data[jj];
      do {
        exitg2 = 0;
        jj++;
        if (jj + 1 > k) {
          exitg2 = 1;
        } else {
          absx = fabs(x / 2.0);
          if (absx <= 2.2250738585072014E-308) {
            absx = 4.94065645841247E-324;
          } else {
            frexp(absx, &exponent);
            absx = ldexp(1.0, exponent - 53);
          }
          if (!(fabs(x - d_data[jj]) < absx)) {
            exitg2 = 1;
          }
        }
      } while (exitg2 == 0);
      ii++;
      d_data[ii] = x;
    }
    if (1 > ii + 1) {
      b_i = 0;
    } else {
      b_i = ii + 1;
    }
    i1 = d->size[0];
    d->size[0] = b_i;
    emxEnsureCapacity_real_T(d, i1);
    d_data = d->data;
    if (b_i == 0) {
      res1->size[0] = arg1->size[0];
      res1->size[1] = 0;
    } else {
      i1 = res1->size[0] * res1->size[1];
      res1->size[0] = (int)fmin(arg1->size[0], arg1->size[1]);
      res1->size[1] = b_i;
      emxEnsureCapacity_real_T(res1, i1);
      res1_data = res1->data;
      ii = (int)fmin(arg1->size[0], arg1->size[1]) * b_i;
      for (i1 = 0; i1 < ii; i1++) {
        res1_data[i1] = 0.0;
      }
      if (arg1->size[0] >= arg1->size[1]) {
        for (k = 0; k < b_i; k++) {
          absx = fmax(1.0, d_data[k] + 1.0);
          x = fmin(n, (double)m + d_data[k]);
          if (x < absx) {
            i->size[0] = 1;
            i->size[1] = 0;
          } else if (absx == absx) {
            i1 = i->size[0] * i->size[1];
            i->size[0] = 1;
            ii = (int)(x - absx);
            i->size[1] = ii + 1;
            emxEnsureCapacity_real_T(i, i1);
            i_data = i->data;
            for (i1 = 0; i1 <= ii; i1++) {
              i_data[i1] = absx + (double)i1;
            }
          } else {
            eml_float_colon(absx, x, i);
            i_data = i->data;
          }
          diagk(arg1, d_data[k], a);
          a_data = a->data;
          ii = a->size[0];
          for (i1 = 0; i1 < ii; i1++) {
            res1_data[((int)i_data[i1] + res1->size[0] * k) - 1] = a_data[i1];
          }
        }
      } else {
        for (k = 0; k < b_i; k++) {
          absx = fmax(1.0, 1.0 - d_data[k]);
          x = fmin(m, (double)n - d_data[k]);
          if (x < absx) {
            i->size[0] = 1;
            i->size[1] = 0;
          } else if (absx == absx) {
            i1 = i->size[0] * i->size[1];
            i->size[0] = 1;
            ii = (int)(x - absx);
            i->size[1] = ii + 1;
            emxEnsureCapacity_real_T(i, i1);
            i_data = i->data;
            for (i1 = 0; i1 <= ii; i1++) {
              i_data[i1] = absx + (double)i1;
            }
          } else {
            eml_float_colon(absx, x, i);
            i_data = i->data;
          }
          diagk(arg1, d_data[k], a);
          a_data = a->data;
          ii = a->size[0];
          for (i1 = 0; i1 < ii; i1++) {
            res1_data[((int)i_data[i1] + res1->size[0] * k) - 1] = a_data[i1];
          }
        }
      }
    }
  }
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&idx);
  emxFree_real_T(&a);
  emxFree_int32_T(&j);
  emxFree_real_T(&i);
  emxFree_real_T(&d);
}

/* End of code generation (spdiags.c) */
