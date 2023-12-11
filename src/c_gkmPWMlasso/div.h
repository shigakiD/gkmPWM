/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * div.h
 *
 * Code generation for function 'div'
 *
 */

#ifndef DIV_H
#define DIV_H

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
void v_binary_expand_op(emxArray_real_T *b, const emxArray_real_T *b_b,
                        const emxArray_real_T *fit);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (div.h) */
