/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 */

#ifndef REPMAT_H
#define REPMAT_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_repmat(const double a[4], double varargin_1, emxArray_real_T *b);

void c_repmat(const emxArray_real_T *a, emxArray_real_T *b);

void repmat(const double a[4], double b[16]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (repmat.h) */
