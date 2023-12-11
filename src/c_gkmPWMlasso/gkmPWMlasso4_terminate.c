/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso4_terminate.c
 *
 * Code generation for function 'gkmPWMlasso4_terminate'
 *
 */

/* Include files */
#include "gkmPWMlasso4_terminate.h"
#include "gkmPWMlasso4_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWMlasso4_terminate(void)
{
  omp_destroy_nest_lock(&gkmPWMlasso4_nestLockGlobal);
  isInitialized_gkmPWMlasso4 = false;
}

/* End of code generation (gkmPWMlasso4_terminate.c) */
