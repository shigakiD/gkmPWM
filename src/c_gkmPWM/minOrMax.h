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
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_maximum(const double x[9], double *ex, int *idx);

creal_T b_minimum(const creal_T x[4]);

void c_maximum(const double x[2], double *ex, int *idx);

void d_maximum(const emxArray_real_T *x, double *ex, int *idx);

double e_maximum(const double x[4]);

void f_maximum(const emxArray_real_T *x, emxArray_real_T *ex);

void g_maximum(const emxArray_real_T *x, emxArray_real_T *ex);

double h_maximum(const emxArray_real_T *x);

void i_maximum(const creal_T x[4], creal_T *ex, int *idx);

double j_maximum(const double x[9]);

double maximum(const emxArray_real_T *x);

double minimum(const double x[2]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (minOrMax.h) */