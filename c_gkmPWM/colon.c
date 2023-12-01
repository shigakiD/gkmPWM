/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * colon.c
 *
 * Code generation for function 'colon'
 *
 */

/* Include files */
#include "colon.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void eml_float_colon(double d, double b, emxArray_real_T *y)
{
  double apnd;
  double cdiff;
  double ndbl;
  double *y_data;
  int k;
  int n;
  int nm1d2;
  ndbl = floor((b - 1.0) / d + 0.5);
  apnd = ndbl * d + 1.0;
  cdiff = apnd - b;
  if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(1.0, fabs(b))) {
    ndbl++;
    apnd = b;
  } else if (cdiff > 0.0) {
    apnd = (ndbl - 1.0) * d + 1.0;
  } else {
    ndbl++;
  }
  if (ndbl >= 0.0) {
    n = (int)ndbl;
  } else {
    n = 0;
  }
  nm1d2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity_real_T(y, nm1d2);
  y_data = y->data;
  if (n > 0) {
    y_data[0] = 1.0;
    if (n > 1) {
      y_data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        ndbl = ((double)k + 1.0) * d;
        y_data[k + 1] = ndbl + 1.0;
        y_data[(n - k) - 2] = apnd - ndbl;
      }
      if (nm1d2 << 1 == n - 1) {
        y_data[nm1d2] = (apnd + 1.0) / 2.0;
      } else {
        ndbl = (double)nm1d2 * d;
        y_data[nm1d2] = ndbl + 1.0;
        y_data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }
}

/* End of code generation (colon.c) */
