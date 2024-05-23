/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso_initialize.c
 *
 * Code generation for function 'gkmPWMlasso_initialize'
 *
 */

/* Include files */
#include "gkmPWMlasso_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include "gkmPWMlasso_data.h"
#include <string.h>

/* Function Definitions */
void gkmPWMlasso_initialize(void)
{
  omp_init_nest_lock(&gkmPWMlasso_nestLockGlobal);
  omp_set_num_threads(1);  
  c_eml_rand_mt19937ar_stateful_i();
  filedata_init();
  isInitialized_gkmPWMlasso = true;
}

/* End of code generation (gkmPWMlasso_initialize.c) */
