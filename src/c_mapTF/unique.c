/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * unique.c
 *
 * Code generation for function 'unique'
 *
 */

/* Include files */
#include "unique.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void unique_vector(const emxArray_real_T *a, emxArray_real_T *b)
{
  emxArray_int32_T *idx;
  emxArray_int32_T *iwork;
  const double *a_data;
  double absx;
  double x;
  double *b_data;
  int b_i;
  int exitg2;
  int exponent;
  int i;
  int i2;
  int j;
  int k;
  int kEnd;
  int n;
  int na;
  int nb;
  int pEnd;
  int q;
  int qEnd;
  int *idx_data;
  int *iwork_data;
  bool exitg1;
  a_data = a->data;
  emxInit_int32_T(&idx, 1);
  na = a->size[0];
  n = a->size[0] + 1;
  i = idx->size[0];
  idx->size[0] = a->size[0];
  emxEnsureCapacity_int32_T(idx, i);
  idx_data = idx->data;
  b_i = a->size[0];
  for (i = 0; i < b_i; i++) {
    idx_data[i] = 0;
  }
  if (a->size[0] != 0) {
    emxInit_int32_T(&iwork, 1);
    i = iwork->size[0];
    iwork->size[0] = a->size[0];
    emxEnsureCapacity_int32_T(iwork, i);
    iwork_data = iwork->data;
    i = a->size[0] - 1;
    for (k = 1; k <= i; k += 2) {
      if ((a_data[k - 1] <= a_data[k]) || rtIsNaN(a_data[k])) {
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
    b_i = 2;
    while (b_i < n - 1) {
      i2 = b_i << 1;
      j = 1;
      for (pEnd = b_i + 1; pEnd < n; pEnd = qEnd + b_i) {
        nb = j;
        q = pEnd - 1;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }
        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          absx = a_data[idx_data[q] - 1];
          i = idx_data[nb - 1];
          if ((a_data[i - 1] <= absx) || rtIsNaN(absx)) {
            iwork_data[k] = i;
            nb++;
            if (nb == pEnd) {
              while (q + 1 < qEnd) {
                k++;
                iwork_data[k] = idx_data[q];
                q++;
              }
            }
          } else {
            iwork_data[k] = idx_data[q];
            q++;
            if (q + 1 == qEnd) {
              while (nb < pEnd) {
                k++;
                iwork_data[k] = idx_data[nb - 1];
                nb++;
              }
            }
          }
          k++;
        }
        for (k = 0; k < kEnd; k++) {
          idx_data[(j + k) - 1] = iwork_data[k];
        }
        j = qEnd;
      }
      b_i = i2;
    }
    emxFree_int32_T(&iwork);
  }
  i = b->size[0];
  b->size[0] = a->size[0];
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  for (k = 0; k < na; k++) {
    b_data[k] = a_data[idx_data[k] - 1];
  }
  emxFree_int32_T(&idx);
  k = 0;
  while ((k + 1 <= na) && rtIsInf(b_data[k]) && (b_data[k] < 0.0)) {
    k++;
  }
  i2 = k;
  k = a->size[0];
  while ((k >= 1) && rtIsNaN(b_data[k - 1])) {
    k--;
  }
  pEnd = a->size[0] - k;
  exitg1 = false;
  while ((!exitg1) && (k >= 1)) {
    absx = b_data[k - 1];
    if (rtIsInf(absx) && (absx > 0.0)) {
      k--;
    } else {
      exitg1 = true;
    }
  }
  b_i = (a->size[0] - k) - pEnd;
  nb = -1;
  if (i2 > 0) {
    nb = 0;
  }
  while (i2 + 1 <= k) {
    x = b_data[i2];
    do {
      exitg2 = 0;
      i2++;
      if (i2 + 1 > k) {
        exitg2 = 1;
      } else {
        absx = fabs(x / 2.0);
        if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
          if (absx <= 2.2250738585072014E-308) {
            absx = 4.94065645841247E-324;
          } else {
            frexp(absx, &exponent);
            absx = ldexp(1.0, exponent - 53);
          }
        } else {
          absx = rtNaN;
        }
        if ((!(fabs(x - b_data[i2]) < absx)) &&
            ((!rtIsInf(b_data[i2])) || (!rtIsInf(x)) ||
             ((b_data[i2] > 0.0) != (x > 0.0)))) {
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0);
    nb++;
    b_data[nb] = x;
  }
  if (b_i > 0) {
    nb++;
    b_data[nb] = b_data[k];
  }
  i2 = k + b_i;
  for (j = 0; j < pEnd; j++) {
    b_data[(nb + j) + 1] = b_data[i2 + j];
  }
  nb += pEnd;
  i = b->size[0];
  if (1 > nb + 1) {
    b->size[0] = 0;
  } else {
    b->size[0] = nb + 1;
  }
  emxEnsureCapacity_real_T(b, i);
}

/* End of code generation (unique.c) */
