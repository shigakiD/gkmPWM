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
    *x = a.re;
    *y = b.re;
}

/* End of code generation (relop.c) */
