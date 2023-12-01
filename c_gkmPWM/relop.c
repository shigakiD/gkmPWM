/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * relop.c
 *
 * Code generation for function 'relop'
 *
 */

/* Include files */
#include "relop.h"
#include "gkmPWM_rtwutil.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void absRelopProxies(const creal_T a, const creal_T b, double *x, double *y)
{
  double absx;
  int b_exponent;
  int c_exponent;
  int exponent;
  bool SCALEA;
  bool SCALEB;
  if ((fabs(a.re) > 8.9884656743115785E+307) ||
      (fabs(a.im) > 8.9884656743115785E+307)) {
    SCALEA = true;
  } else {
    SCALEA = false;
  }
  if ((fabs(b.re) > 8.9884656743115785E+307) ||
      (fabs(b.im) > 8.9884656743115785E+307)) {
    SCALEB = true;
  } else {
    SCALEB = false;
  }
  if (SCALEA || SCALEB) {
    *x = rt_hypotd(a.re / 2.0, a.im / 2.0);
    *y = rt_hypotd(b.re / 2.0, b.im / 2.0);
  } else {
    *x = rt_hypotd(a.re, a.im);
    *y = rt_hypotd(b.re, b.im);
  }
  absx = *y / 2.0;
  if (absx <= 2.2250738585072014E-308) {
    absx = 4.94065645841247E-324;
  } else {
    frexp(absx, &exponent);
    absx = ldexp(1.0, exponent - 53);
  }
  if (fabs(*y - *x) < absx) {
    *x = atan2(a.im, a.re);
    *y = atan2(b.im, b.re);
    absx = fabs(*y / 2.0);
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &b_exponent);
      absx = ldexp(1.0, b_exponent - 53);
    }
    if (fabs(*y - *x) < absx) {
      if (a.re != b.re) {
        if (*x >= 0.0) {
          *x = b.re;
          *y = a.re;
        } else {
          *x = a.re;
          *y = b.re;
        }
      } else if (a.re < 0.0) {
        *x = b.im;
        *y = a.im;
      } else {
        *x = a.im;
        *y = b.im;
      }
      absx = fabs(*y / 2.0);
      if (absx <= 2.2250738585072014E-308) {
        absx = 4.94065645841247E-324;
      } else {
        frexp(absx, &c_exponent);
        absx = ldexp(1.0, c_exponent - 53);
      }
      if (fabs(*y - *x) < absx) {
        *x = 0.0;
        *y = 0.0;
      }
    }
  }
}

/* End of code generation (relop.c) */
