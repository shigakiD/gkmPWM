/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PWM2kmers_norc.h
 *
 * Code generation for function 'PWM2kmers_norc'
 *
 */

#ifndef PWM2KMERS_NORC_H
#define PWM2KMERS_NORC_H

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
void PWM2kmers_norc(const emxArray_real_T *mat, const double negmat[16],
                    const emxArray_real_T *c, const emxArray_real_T *s,
                    const emxArray_real_T *ind, const emxArray_real_T *indloc,
                    const emxArray_real_T *x, double l, double k,
                    emxArray_real_T *kweig);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (PWM2kmers_norc.h) */
