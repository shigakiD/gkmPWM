/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * conncomp.h
 *
 * Code generation for function 'conncomp'
 *
 */

#ifndef CONNCOMP_H
#define CONNCOMP_H

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
void graph_conncomp(const emxArray_int32_T *G_Underlying_Ir,
                    const emxArray_int32_T *G_Underlying_Jc,
                    emxArray_real_T *bins);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (conncomp.h) */
