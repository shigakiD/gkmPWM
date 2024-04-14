/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM_data.h
 *
 * Code generation for function 'gkmPWM_data'
 *
 */

#ifndef GKMPWM_DATA_H
#define GKMPWM_DATA_H

/* Include files */
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

/* Variable Declarations */
extern unsigned int state[625];
extern double freq;
extern bool freq_not_empty;
extern omp_nest_lock_t gkmPWM_nestLockGlobal;
extern const bool bv[128];
extern bool isInitialized_gkmPWM;

#endif
/* End of code generation (gkmPWM_data.h) */
