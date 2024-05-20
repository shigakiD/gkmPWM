/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmweights_initialize.c
 *
 * Code generation for function 'getgkmweights_initialize'
 *
 */

/* Include files */
#include "getgkmweights_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include "getgkmweights_data.h"

/* Function Definitions */
void getgkmweights_initialize(void)
{
  omp_init_nest_lock(&getgkmweights_nestLockGlobal);
  omp_set_num_threads(1);
  c_eml_rand_mt19937ar_stateful_i();
  filedata_init();
  isInitialized_getgkmweights = true;
}

/* End of code generation (getgkmweights_initialize.c) */
