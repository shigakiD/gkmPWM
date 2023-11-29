/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * std.c
 *
 * Code generation for function 'std'
 *
 */

/* Include files */
#include "std.h"
#include "blockedSummation.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include "cblas.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
double b_std(const emxArray_real_T *x)
{
  emxArray_real_T *absdiff;
  const double *x_data;
  double xbar;
  double y;
  double *absdiff_data;
  int k;
  int n;
  x_data = x->data;
  n = x->size[0];
  if (x->size[0] == 0) {
    y = 0.0;
  } else if (x->size[0] == 1) {
    y = 0.0;
  } else {
    emxInit_real_T(&absdiff, 1);
    xbar = blockedSummation(x, x->size[0]) / (double)x->size[0];
    k = absdiff->size[0];
    absdiff->size[0] = x->size[0];
    emxEnsureCapacity_real_T(absdiff, k);
    absdiff_data = absdiff->data;
    for (k = 0; k < n; k++) {
      absdiff_data[k] = fabs(x_data[k] - xbar);
    }
    y = cblas_dnrm2((blasint)x->size[0], &absdiff_data[0], (blasint)1);
    y /= sqrt((double)x->size[0] - 1.0);
    emxFree_real_T(&absdiff);
  }
  return y;
}

/* End of code generation (std.c) */
