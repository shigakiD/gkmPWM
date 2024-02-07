/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nchoosek.c
 *
 * Code generation for function 'nchoosek'
 *
 */

/* Include files */
#include "nchoosek.h"
#include "getgkmweights_emxutil.h"
#include "getgkmweights_rtwutil.h"
#include "getgkmweights_types.h"
#include <math.h>

/* Function Definitions */
/*
 *
 */
double b_nchoosek(double x, double k)
{
  double b_k;
  double nmk;
  double y;
  int i;
  int j;
  if (k <= x) {
    b_k = k;
    if (k > x / 2.0) {
      b_k = x - k;
    }
    if (b_k > 1000.0) {
      y = 1.7976931348623157E+308;
    } else {
      y = 1.0;
      nmk = x - b_k;
      i = (int)b_k;
      for (j = 0; j < i; j++) {
        y *= (((double)j + 1.0) + nmk) / ((double)j + 1.0);
      }
      y = rt_roundd(y);
    }
  }
  return y;
}

/*
 *
 */
void nchoosek(const emxArray_real_T *x, double k, emxArray_real_T *y)
{
  emxArray_int32_T *comb;
  const double *x_data;
  double b_k;
  double b_y;
  double nmk;
  double *y_data;
  int combj;
  int kint;
  int n;
  int nmkpi;
  int nrows_tmp;
  int row;
  int yk;
  int *comb_data;
  x_data = x->data;
  if (x->size[1] == 1) {
    if (k <= x_data[0]) {
      b_k = k;
      if (k > x_data[0] / 2.0) {
        b_k = x_data[0] - k;
      }
      if (b_k > 1000.0) {
        b_y = 1.7976931348623157E+308;
      } else {
        b_y = 1.0;
        nmk = x_data[0] - b_k;
        yk = (int)b_k;
        for (n = 0; n < yk; n++) {
          b_y *= (((double)n + 1.0) + nmk) / ((double)n + 1.0);
        }
        b_y = rt_roundd(b_y);
      }
      yk = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, yk);
      y_data = y->data;
      y_data[0] = b_y;
    }
  } else if (k > x->size[1]) {
    y->size[0] = 0;
    y->size[1] = (int)k;
  } else {
    kint = (int)floor(k);
    b_k = k;
    if (k > (double)x->size[1] / 2.0) {
      b_k = (double)x->size[1] - k;
    }
    if (b_k > 1000.0) {
      b_y = 1.7976931348623157E+308;
    } else {
      b_y = 1.0;
      nmk = (double)x->size[1] - b_k;
      yk = (int)b_k;
      for (n = 0; n < yk; n++) {
        b_y *= (((double)n + 1.0) + nmk) / ((double)n + 1.0);
      }
      b_y = rt_roundd(b_y);
    }
    emxInit_int32_T(&comb, 2);
    nrows_tmp = (int)b_y;
    yk = y->size[0] * y->size[1];
    y->size[0] = (int)b_y;
    y->size[1] = kint;
    emxEnsureCapacity_real_T(y, yk);
    y_data = y->data;
    if (kint < 1) {
      n = 0;
    } else {
      n = kint;
    }
    yk = comb->size[0] * comb->size[1];
    comb->size[0] = 1;
    comb->size[1] = n;
    emxEnsureCapacity_int32_T(comb, yk);
    comb_data = comb->data;
    if (n > 0) {
      comb_data[0] = 1;
      yk = 1;
      for (nmkpi = 2; nmkpi <= n; nmkpi++) {
        yk++;
        comb_data[nmkpi - 1] = yk;
      }
    }
    n = kint - 1;
    nmkpi = x->size[1];
    for (row = 0; row < nrows_tmp; row++) {
      for (yk = 0; yk < kint; yk++) {
        y_data[row + y->size[0] * yk] = x_data[comb_data[yk] - 1];
      }
      if (n + 1 > 0) {
        yk = comb_data[n];
        combj = comb_data[n] + 1;
        comb_data[n]++;
        if (yk + 1 < nmkpi) {
          yk = n + 2;
          for (n = yk; n <= kint; n++) {
            combj++;
            comb_data[n - 1] = combj;
          }
          n = kint - 1;
          nmkpi = x->size[1];
        } else {
          n--;
          nmkpi--;
        }
      }
    }
    emxFree_int32_T(&comb);
  }
}

/* End of code generation (nchoosek.c) */
