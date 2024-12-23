/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nullAssignment.h
 *
 * Code generation for function 'nullAssignment'
 *
 */

#ifndef NULLASSIGNMENT_H
#define NULLASSIGNMENT_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_nullAssignment(emxArray_real_T *x, int idx);

void c_nullAssignment(emxArray_real_T *x, const emxArray_int32_T *idx);

void d_nullAssignment(emxArray_real_T *x, int idx);

void e_nullAssignment(emxArray_real_T *x, const emxArray_int32_T *idx);

void nullAssignment(emxArray_real_T *x, const emxArray_int32_T *idx);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (nullAssignment.h) */
