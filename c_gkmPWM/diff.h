/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff.h
 *
 * Code generation for function 'diff'
 *
 */

#ifndef DIFF_H
#define DIFF_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_diff(const double x[5], double y[4]);

void c_diff(const double x[10], double y[9]);

void d_diff(const emxArray_boolean_T *x, emxArray_real_T *y);

void diff(const emxArray_real_T *x, emxArray_real_T *y);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (diff.h) */
