/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmcounts.h
 *
 * Code generation for function 'getgkmcounts'
 *
 */

#ifndef GETGKMCOUNTS_H
#define GETGKMCOUNTS_H

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
void getgkmcounts(const emxArray_char_T *filename, double l, double k,
                  double lk, bool RC, emxArray_real_T *gkmc, double *GCpos,
                  double *GCneg, double mat[16], double mat2[16]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (getgkmcounts.h) */