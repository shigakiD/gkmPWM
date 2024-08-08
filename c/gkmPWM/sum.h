/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.h
 *
 * Code generation for function 'sum'
 *
 */

#ifndef SUM_H
#define SUM_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_sum(const emxArray_real_T *x, emxArray_real_T *y);

void c_sum(const emxArray_real_T *x, emxArray_real_T *y);

void d_sum(const double x[16], double y[4]);

void e_sum(const creal_T x[16], creal_T y[4]);

void f_sum(const creal_T x_data[], const int x_size[2], creal_T y_data[],
           int y_size[2]);

double sum(const emxArray_real_T *x);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (sum.h) */
