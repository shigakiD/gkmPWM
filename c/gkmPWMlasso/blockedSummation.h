/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * blockedSummation.h
 *
 * Code generation for function 'blockedSummation'
 *
 */

#ifndef BLOCKEDSUMMATION_H
#define BLOCKEDSUMMATION_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
double blockedSummation(const emxArray_real_T *x, int vlen);

void colMajorFlatIter(const emxArray_real_T *x, int vlen, emxArray_real_T *y);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (blockedSummation.h) */
