/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xnrm2.c
 *
 * Code generation for function 'xnrm2'
 *
 */

/* Include files */
#include "xnrm2.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
double b_xnrm2(int n, const double x[3])
{
  double absxk;
  double scale;
  double t;
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[1]);
    } else {
      scale = 3.3121686421112381E-170;
      for (k = 0; k < 2; k++) {
        absxk = fabs(x[k + 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

/*
 *
 */
double xnrm2(int n, const double x[16], int ix0)
{
  double absxk;
  double scale;
  double t;
  double y;
  int k;
  int kend;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = ix0 + 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

/* End of code generation (xnrm2.c) */
