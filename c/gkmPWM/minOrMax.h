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
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_maximum(const double x[9], double *ex, int *idx);

double b_minimum(const double x[2]);

void c_maximum(const double x[2], double *ex, int *idx);

creal_T c_minimum(const creal_T x[4]);

void d_maximum(const emxArray_real_T *x, double *ex, int *idx);

double e_maximum(const double x[4]);

void f_maximum(const emxArray_real_T *x, emxArray_real_T *ex);

void g_maximum(const emxArray_real_T *x, emxArray_real_T *ex);

double h_maximum(const emxArray_real_T *x);

double i_maximum(const double x[9]);

double maximum(const emxArray_real_T *x);

double minimum(const emxArray_real_T *x);

void q_binary_expand_op(emxArray_real_T *Rd, const emxArray_real_T *CT,
                        const emxArray_real_T *ct);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (minOrMax.h) */
