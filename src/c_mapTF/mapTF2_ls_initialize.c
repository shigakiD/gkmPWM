/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF2_ls_initialize.c
 *
 * Code generation for function 'mapTF2_ls_initialize'
 *
 */

/* Include files */
#include "mapTF2_ls_initialize.h"
#include "CoderTimeAPI.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include "mapTF2_ls_data.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"
#include <string.h>

/* Function Definitions */
void mapTF2_ls_initialize(void)
{
  omp_init_nest_lock(&mapTF2_ls_nestLockGlobal);
  savedTime_not_empty_init();
  freq_not_empty_init();
  filedata_init();
  c_eml_rand_mt19937ar_stateful_i();
  isInitialized_mapTF2_ls = true;
}

/* End of code generation (mapTF2_ls_initialize.c) */
