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
#include "mapTF2_ls_emxutil.h"
#include "mapTF2_ls_rtwutil.h"
#include "mapTF2_ls_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static double nCk(double n, double k);

/* Function Definitions */
/*
 *
 */
static double nCk(double n, double k)
{
  double nmk;
  double y;
  int i;
  int j;
  if ((!rtIsInf(n)) && (!rtIsNaN(n)) && ((!rtIsInf(k)) && (!rtIsNaN(k)))) {
    if (k > n / 2.0) {
      k = n - k;
    }
    if (k > 1000.0) {
      y = rtInf;
    } else {
      y = 1.0;
      nmk = n - k;
      i = (int)k;
      for (j = 0; j < i; j++) {
        y *= (((double)j + 1.0) + nmk) / ((double)j + 1.0);
      }
      y = rt_roundd_snf(y);
    }
  } else {
    y = rtNaN;
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
  double *y_data;
  int combj;
  int kint;
  int n;
  int nmkpi;
  int nrows;
  int row;
  int yk;
  int *comb_data;
  x_data = x->data;
  if (x->size[1] == 1) {
    if (!(k > x_data[0])) {
      yk = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, yk);
      y_data = y->data;
      y_data[0] = nCk(x_data[0], k);
    }
  } else if (k > x->size[1]) {
    y->size[0] = 0;
    y->size[1] = (int)k;
  } else {
    emxInit_int32_T(&comb, 2);
    kint = (int)floor(k);
    nrows = (int)floor(nCk(x->size[1], k));
    yk = y->size[0] * y->size[1];
    y->size[0] = nrows;
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
    for (row = 0; row < nrows; row++) {
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
