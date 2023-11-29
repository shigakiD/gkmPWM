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
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
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
  int exitg1;
  int exponent;
  int i;
  int i1;
  int i2;
  int j;
  int k;
  int kEnd;
  int n;
  int na;
  int p;
  int pEnd;
  int q;
  int qEnd;
  int *idx_data;
  int *iwork_data;
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
    b_i = 2;
    while (b_i < n - 1) {
      i2 = b_i << 1;
      j = 1;
      for (pEnd = b_i + 1; pEnd < n; pEnd = qEnd + b_i) {
        p = j;
        q = pEnd;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }
        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          i = idx_data[q - 1];
          i1 = idx_data[p - 1];
          if (a_data[i1 - 1] <= a_data[i - 1]) {
            iwork_data[k] = i1;
            p++;
            if (p == pEnd) {
              while (q < qEnd) {
                k++;
                iwork_data[k] = idx_data[q - 1];
                q++;
              }
            }
          } else {
            iwork_data[k] = i;
            q++;
            if (q == qEnd) {
              while (p < pEnd) {
                k++;
                iwork_data[k] = idx_data[p - 1];
                p++;
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
  b_i = 0;
  k = 1;
  while (k <= na) {
    x = b_data[k - 1];
    do {
      exitg1 = 0;
      k++;
      if (k > na) {
        exitg1 = 1;
      } else {
        absx = fabs(x / 2.0);
        if (absx <= 2.2250738585072014E-308) {
          absx = 4.94065645841247E-324;
        } else {
          frexp(absx, &exponent);
          absx = ldexp(1.0, exponent - 53);
        }
        if (fabs(x - b_data[k - 1]) >= absx) {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
    b_i++;
    b_data[b_i - 1] = x;
  }
  i = b->size[0];
  if (1 > b_i) {
    b->size[0] = 0;
  } else {
    b->size[0] = b_i;
  }
  emxEnsureCapacity_real_T(b, i);
}

/* End of code generation (unique.c) */
