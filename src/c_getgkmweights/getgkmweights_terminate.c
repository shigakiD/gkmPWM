/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmweights_terminate.c
 *
 * Code generation for function 'getgkmweights_terminate'
 *
 */

/* Include files */
#include "getgkmweights_terminate.h"
#include "getgkmweights_data.h"

/* Function Definitions */
void getgkmweights_terminate(void)
{
  omp_destroy_nest_lock(&getgkmweights_nestLockGlobal);
  isInitialized_getgkmweights = false;
}

/* End of code generation (getgkmweights_terminate.c) */
