/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso_data.h
 *
 * Code generation for function 'gkmPWMlasso_data'
 *
 */

#ifndef GKMPWMLASSO4_DATA_H
#define GKMPWMLASSO4_DATA_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

/* Variable Declarations */
extern unsigned int state[625];
extern omp_nest_lock_t gkmPWMlasso_nestLockGlobal;
extern const bool bv[128];
extern const struct_T r;
extern bool isInitialized_gkmPWMlasso;

#endif
/* End of code generation (gkmPWMlasso_data.h) */