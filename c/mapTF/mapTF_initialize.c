/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF_initialize.c
 *
 * Code generation for function 'mapTF_initialize'
 *
 */

/* Include files */
#include "mapTF_initialize.h"
#include "CoderTimeAPI.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "fileManager.h"
#include "mapTF_data.h"
#include "rt_nonfinite.h"
#include "timeKeeper.h"
#include <string.h>

/* Function Definitions */
void mapTF_initialize(void)
{
  savedTime_not_empty_init();
  freq_not_empty_init();
  filedata_init();
  c_eml_rand_mt19937ar_stateful_i();
  isInitialized_mapTF = true;
}

/* End of code generation (mapTF_initialize.c) */
