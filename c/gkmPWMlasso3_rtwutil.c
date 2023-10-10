/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso3_rtwutil.c
 *
 * Code generation for function 'gkmPWMlasso3_rtwutil'
 *
 */

/* Include files */
#include "gkmPWMlasso3_rtwutil.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
double rt_roundd(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }
  return y;
}

/* End of code generation (gkmPWMlasso3_rtwutil.c) */
