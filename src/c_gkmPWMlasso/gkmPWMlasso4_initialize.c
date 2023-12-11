/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso4_initialize.c
 *
 * Code generation for function 'gkmPWMlasso4_initialize'
 *
 */

/* Include files */
#include "gkmPWMlasso4_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include "gkmPWMlasso4_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWMlasso4_initialize(void)
{
  omp_init_nest_lock(&gkmPWMlasso4_nestLockGlobal);
  c_eml_rand_mt19937ar_stateful_i();
  filedata_init();
  isInitialized_gkmPWMlasso4 = true;
}

/* End of code generation (gkmPWMlasso4_initialize.c) */
