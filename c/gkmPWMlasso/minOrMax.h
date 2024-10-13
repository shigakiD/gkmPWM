/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minOrMax.h
 *
 * Code generation for function 'minOrMax'
 *
 */

#ifndef MINORMAX_H
#define MINORMAX_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
double b_maximum(const double x[2]);

void c_maximum(const emxArray_real_T *x, double *ex, int *idx);

void d_maximum(const emxArray_real_T *x, emxArray_real_T *ex);

void e_maximum(const emxArray_real_T *x, double *ex, int *idx);

double maximum(const emxArray_real_T *x);

void minimum(const emxArray_real_T *x, emxArray_real_T *ex);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (minOrMax.h) */
