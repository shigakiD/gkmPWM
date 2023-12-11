/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF2_ls_terminate.c
 *
 * Code generation for function 'mapTF2_ls_terminate'
 *
 */

/* Include files */
#include "mapTF2_ls_terminate.h"
#include "mapTF2_ls_data.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void mapTF2_ls_terminate(void)
{
  omp_destroy_nest_lock(&mapTF2_ls_nestLockGlobal);
  isInitialized_mapTF2_ls = false;
}

/* End of code generation (mapTF2_ls_terminate.c) */
