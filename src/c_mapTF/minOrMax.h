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
#include "mapTF_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_maximum(const emxArray_real_T *x, emxArray_real_T *ex);

void c_maximum(const emxArray_real_T *x, double *ex, int *idx);

int d_binary_expand_op(emxArray_real_T *mat, int i1, int i,
                       const emxArray_real_T *TRANS, int i2,
                       const emxArray_real_T *C, int j,
                       const emxArray_real_T *pwm_prob);

int e_binary_expand_op(emxArray_real_T *mat, int i1, int i,
                       const emxArray_real_T *TRANS, int i2,
                       const emxArray_real_T *neg, double y);

double maximum(const emxArray_real_T *x);

void minimum(const emxArray_real_T *x, emxArray_real_T *ex);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (minOrMax.h) */
