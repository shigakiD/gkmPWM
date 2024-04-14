/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM_initialize.c
 *
 * Code generation for function 'gkmPWM_initialize'
 *
 */

/* Include files */
#include "gkmPWM_initialize.h"
#include "CoderTimeAPI.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include "gkmPWM_data.h"
#include "timeKeeper.h"
#include <string.h>

/* Function Definitions */
void gkmPWM_initialize(void)
{
  omp_init_nest_lock(&gkmPWM_nestLockGlobal);
  savedTime_not_empty_init();
  freq_not_empty_init();
  c_eml_rand_mt19937ar_stateful_i();
  filedata_init();
  isInitialized_gkmPWM = true;
}

/* End of code generation (gkmPWM_initialize.c) */
