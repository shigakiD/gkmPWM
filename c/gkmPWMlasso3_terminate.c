/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso3_terminate.c
 *
 * Code generation for function 'gkmPWMlasso3_terminate'
 *
 */

/* Include files */
#include "gkmPWMlasso3_terminate.h"
#include "gkmPWMlasso3_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWMlasso3_terminate(void)
{
  omp_destroy_nest_lock(&gkmPWMlasso3_nestLockGlobal);
  isInitialized_gkmPWMlasso3 = false;
}

/* End of code generation (gkmPWMlasso3_terminate.c) */
