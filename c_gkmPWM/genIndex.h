/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * genIndex.h
 *
 * Code generation for function 'genIndex'
 *
 */

#ifndef GENINDEX_H
#define GENINDEX_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_genIndex(double l, double k, emxArray_real_T *c);

void genIndex(double l, double k, emxArray_real_T *c, emxArray_real_T *C,
              emxArray_real_T *b_I, emxArray_real_T *ind, emxArray_real_T *mat,
              double *rcnum);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (genIndex.h) */