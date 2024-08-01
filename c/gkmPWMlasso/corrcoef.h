/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * corrcoef.h
 *
 * Code generation for function 'corrcoef'
 *
 */

#ifndef CORRCOEF_H
#define CORRCOEF_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_corrcoef(const emxArray_real_T *x, emxArray_real_T *b_r);

void corrcoef(const emxArray_real_T *x, const emxArray_real_T *varargin_1,
              double b_r[4]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (corrcoef.h) */
