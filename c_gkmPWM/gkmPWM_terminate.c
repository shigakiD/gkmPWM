/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM_terminate.c
 *
 * Code generation for function 'gkmPWM_terminate'
 *
 */

/* Include files */
#include "gkmPWM_terminate.h"
#include "gkmPWM_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWM_terminate(void)
{
  omp_destroy_nest_lock(&gkmPWM_nestLockGlobal);
  isInitialized_gkmPWM = false;
}

/* End of code generation (gkmPWM_terminate.c) */
