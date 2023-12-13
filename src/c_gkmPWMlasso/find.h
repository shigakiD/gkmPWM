/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find.h
 *
 * Code generation for function 'find'
 *
 */

#ifndef FIND_H
#define FIND_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_binary_expand_op(emxArray_int32_T *idx, const emxArray_real_T *indc,
                        const emxArray_real_T *BY);

void b_eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i);

void e_binary_expand_op(emxArray_int32_T *idx, const emxArray_real_T *loc,
                        const emxArray_real_T *Z, double varargin_4);

void eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (find.h) */
