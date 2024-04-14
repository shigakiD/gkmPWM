/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sqrt.c
 *
 * Code generation for function 'sqrt'
 *
 */

/* Include files */
#include "sqrt.h"
#include "gkmPWM_rtwutil.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_sqrt(creal_T *x)
{
  double absxi;
  double absxr;
  double xi;
  double xr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxr = 0.0;
      absxi = sqrt(-xr);
    } else {
      absxr = sqrt(xr);
      absxi = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxr = sqrt(-xi / 2.0);
      absxi = -absxr;
    } else {
      absxr = sqrt(xi / 2.0);
      absxi = absxr;
    }
  } else {
    absxr = fabs(xr);
    absxi = fabs(xi);
    if ((absxr > 4.4942328371557893E+307) ||
        (absxi > 4.4942328371557893E+307)) {
      absxr *= 0.5;
      absxi = rt_hypotd(absxr, absxi * 0.5);
      if (absxi > absxr) {
        absxr = sqrt(absxi) * sqrt(absxr / absxi + 1.0);
      } else {
        absxr = sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxr = sqrt((rt_hypotd(absxr, absxi) + absxr) * 0.5);
    }
    if (xr > 0.0) {
      absxi = 0.5 * (xi / absxr);
    } else {
      if (xi < 0.0) {
        absxi = -absxr;
      } else {
        absxi = absxr;
      }
      absxr = 0.5 * (xi / absxi);
    }
  }
  x->re = absxr;
  x->im = absxi;
}

/* End of code generation (sqrt.c) */
