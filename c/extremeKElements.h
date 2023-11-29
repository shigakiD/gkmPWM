/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * extremeKElements.h
 *
 * Code generation for function 'extremeKElements'
 *
 */

#ifndef EXTREMEKELEMENTS_H
#define EXTREMEKELEMENTS_H

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
void exkib(const emxArray_real_T *a, int k, emxArray_int32_T *idx,
           emxArray_real_T *b);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (extremeKElements.h) */
