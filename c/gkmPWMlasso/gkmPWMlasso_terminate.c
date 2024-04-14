/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso_terminate.c
 *
 * Code generation for function 'gkmPWMlasso_terminate'
 *
 */

/* Include files */
#include "gkmPWMlasso_terminate.h"
#include "gkmPWMlasso_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWMlasso_terminate(void)
{
  omp_destroy_nest_lock(&gkmPWMlasso_nestLockGlobal);
  isInitialized_gkmPWMlasso = false;
}

/* End of code generation (gkmPWMlasso_terminate.c) */
