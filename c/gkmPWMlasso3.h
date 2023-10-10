/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso3.h
 *
 * Code generation for function 'gkmPWMlasso3'
 *
 */

#ifndef GKMPWMLASSO3_H
#define GKMPWMLASSO3_H

/* Include files */
#include "gkmPWMlasso3_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void gkmPWMlasso3(const emxArray_char_T *varargin_1,
                         const emxArray_char_T *varargin_2, double varargin_3,
                         double varargin_4, double varargin_5,
                         double varargin_6, double varargin_7,
                         double varargin_8, bool varargin_9);

void minus(emxArray_real_T *normvec, const emxArray_real_T *cfile2);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (gkmPWMlasso3.h) */
