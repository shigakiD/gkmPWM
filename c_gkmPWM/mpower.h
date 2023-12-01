/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mpower.h
 *
 * Code generation for function 'mpower'
 *
 */

#ifndef MPOWER_H
#define MPOWER_H

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
void b_mpower(const creal_T a[16], creal_T c[16]);

void c_mpower(const double a[16], double c[16]);

void d_mpower(const creal_T a[9], creal_T c[9]);

void e_mpower(const creal_T a[4], creal_T c[4]);

void mpower(const emxArray_real_T *a, emxArray_real_T *c);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (mpower.h) */
