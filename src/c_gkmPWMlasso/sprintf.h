/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sprintf.h
 *
 * Code generation for function 'sprintf'
 *
 */

#ifndef SPRINTF_H
#define SPRINTF_H

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
void b_sprintf(const emxArray_char_T *varargin_1, int varargin_2,
               int varargin_3, emxArray_char_T *str_Value);

void c_sprintf(const emxArray_char_T *varargin_1, int varargin_2,
               int varargin_3, emxArray_char_T *str_Value);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (sprintf.h) */
