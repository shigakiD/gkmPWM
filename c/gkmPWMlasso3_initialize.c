/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso3_initialize.c
 *
 * Code generation for function 'gkmPWMlasso3_initialize'
 *
 */

/* Include files */
#include "gkmPWMlasso3_initialize.h"
#include "fileManager.h"
#include "gkmPWMlasso3_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWMlasso3_initialize(void)
{
  omp_init_nest_lock(&gkmPWMlasso3_nestLockGlobal);
  filedata_init();
  isInitialized_gkmPWMlasso3 = true;
}

/* End of code generation (gkmPWMlasso3_initialize.c) */
