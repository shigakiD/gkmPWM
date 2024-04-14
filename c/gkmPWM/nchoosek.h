/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nchoosek.h
 *
 * Code generation for function 'nchoosek'
 *
 */

#ifndef NCHOOSEK_H
#define NCHOOSEK_H

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
double b_nchoosek(double x, double k);

void nchoosek(const emxArray_real_T *x, double k, emxArray_real_T *y);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (nchoosek.h) */
