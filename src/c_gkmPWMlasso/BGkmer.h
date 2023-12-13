/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * BGkmer.h
 *
 * Code generation for function 'BGkmer'
 *
 */

#ifndef BGKMER_H
#define BGKMER_H

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
void BGkmer(const double mat[16], double GC, const emxArray_real_T *c,
            double rcnum, double l, double k, bool RC,
            emxArray_real_T *negweights);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (BGkmer.h) */
