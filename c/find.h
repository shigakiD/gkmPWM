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
#include "gkmPWMlasso4_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_binary_expand_op(emxArray_int32_T *idx, const emxArray_real_T *negvec,
                        const emxArray_real_T *BY);

void b_eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i);

void c_eml_find(const bool x[209], int i_data[], int i_size[2]);

void d_eml_find(const bool x[210], int i_data[], int *i_size);

void e_binary_expand_op(emxArray_int32_T *idx, const emxArray_real_T *loc,
                        const emxArray_real_T *Z, double varargin_4);

void eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (find.h) */