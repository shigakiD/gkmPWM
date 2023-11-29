/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * extremeKElements.c
 *
 * Code generation for function 'extremeKElements'
 *
 */

/* Include files */
#include "extremeKElements.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include "sortedInsertion.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void exkib(const emxArray_real_T *a, int k, emxArray_int32_T *idx,
           emxArray_real_T *b)
{
  emxArray_int32_T *itmp;
  emxArray_int32_T *iwork;
  const double *a_data;
  double *b_data;
  int b_i;
  int b_k;
  int i;
  int i1;
  int i2;
  int j;
  int kEnd;
  int n;
  int p;
  int pEnd;
  int q;
  int qEnd;
  int *idx_data;
  int *itmp_data;
  a_data = a->data;
  n = a->size[0];
  i = idx->size[0];
  idx->size[0] = k;
  emxEnsureCapacity_int32_T(idx, i);
  idx_data = idx->data;
  for (i = 0; i < k; i++) {
    idx_data[i] = 0;
  }
  i = b->size[0];
  b->size[0] = k;
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  for (i = 0; i < k; i++) {
    b_data[i] = 0.0;
  }
  if (k != 0) {
    emxInit_int32_T(&itmp, 1);
    emxInit_int32_T(&iwork, 1);
    if ((k > 64) && (k > (a->size[0] >> 6))) {
      n = a->size[0] + 1;
      i = itmp->size[0];
      itmp->size[0] = a->size[0];
      emxEnsureCapacity_int32_T(itmp, i);
      itmp_data = itmp->data;
      b_i = a->size[0];
      for (i = 0; i < b_i; i++) {
        itmp_data[i] = 0;
      }
      if (a->size[0] != 0) {
        i = iwork->size[0];
        iwork->size[0] = a->size[0];
        emxEnsureCapacity_int32_T(iwork, i);
        idx_data = iwork->data;
        i = a->size[0] - 1;
        for (b_k = 1; b_k <= i; b_k += 2) {
          if (a_data[b_k - 1] >= a_data[b_k]) {
            itmp_data[b_k - 1] = b_k;
            itmp_data[b_k] = b_k + 1;
          } else {
            itmp_data[b_k - 1] = b_k + 1;
            itmp_data[b_k] = b_k;
          }
        }
        if ((a->size[0] & 1) != 0) {
          itmp_data[a->size[0] - 1] = a->size[0];
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
            b_k = 0;
            kEnd = qEnd - j;
            while (b_k + 1 <= kEnd) {
              i = itmp_data[q - 1];
              i1 = itmp_data[p - 1];
              if (a_data[i1 - 1] >= a_data[i - 1]) {
                idx_data[b_k] = i1;
                p++;
                if (p == pEnd) {
                  while (q < qEnd) {
                    b_k++;
                    idx_data[b_k] = itmp_data[q - 1];
                    q++;
                  }
                }
              } else {
                idx_data[b_k] = i;
                q++;
                if (q == qEnd) {
                  while (p < pEnd) {
                    b_k++;
                    idx_data[b_k] = itmp_data[p - 1];
                    p++;
                  }
                }
              }
              b_k++;
            }
            for (b_k = 0; b_k < kEnd; b_k++) {
              itmp_data[(j + b_k) - 1] = idx_data[b_k];
            }
            j = qEnd;
          }
          b_i = i2;
        }
      }
      i = idx->size[0];
      idx->size[0] = k;
      emxEnsureCapacity_int32_T(idx, i);
      idx_data = idx->data;
      i = b->size[0];
      b->size[0] = k;
      emxEnsureCapacity_real_T(b, i);
      b_data = b->data;
      for (j = 0; j < k; j++) {
        idx_data[j] = itmp_data[j];
        b_data[j] = a_data[itmp_data[j] - 1];
      }
    } else {
      for (j = 0; j < k; j++) {
        b_i = j;
        sortedInsertion(a_data[j], j + 1, b, &b_i, k, idx);
      }
      i = k + 1;
      for (j = i; j <= n; j++) {
        b_i = k;
        sortedInsertion(a_data[j - 1], j, b, &b_i, k, idx);
      }
    }
    emxFree_int32_T(&iwork);
    emxFree_int32_T(&itmp);
  }
}

/* End of code generation (extremeKElements.c) */
