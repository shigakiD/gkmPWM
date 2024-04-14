/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF_terminate.c
 *
 * Code generation for function 'mapTF_terminate'
 *
 */

/* Include files */
#include "mapTF_terminate.h"
#include "mapTF_data.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void mapTF_terminate(void)
{
  omp_destroy_nest_lock(&mapTF_nestLockGlobal);
  isInitialized_mapTF = false;
}

/* End of code generation (mapTF_terminate.c) */
