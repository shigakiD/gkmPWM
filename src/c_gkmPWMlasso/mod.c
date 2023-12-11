/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mod.c
 *
 * Code generation for function 'mod'
 *
 */

/* Include files */
#include "mod.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
double b_mod(double x, double y)
{
  double b_r;
  double q;
  bool rEQ0;
  b_r = x;
  if (y == 0.0) {
    if (x == 0.0) {
      b_r = y;
    }
  } else if (x == 0.0) {
    b_r = 0.0 / y;
  } else {
    b_r = fmod(x, y);
    rEQ0 = (b_r == 0.0);
    if ((!rEQ0) && (y > floor(y))) {
      q = fabs(x / y);
      rEQ0 = (fabs(q - floor(q + 0.5)) <= 2.2204460492503131E-16 * q);
    }
    if (rEQ0) {
      b_r = 0.0;
    } else if ((x < 0.0) != (y < 0.0)) {
      b_r += y;
    }
  }
  return b_r;
}

/* End of code generation (mod.c) */
