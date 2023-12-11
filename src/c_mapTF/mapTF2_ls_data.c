/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF2_ls_data.c
 *
 * Code generation for function 'mapTF2_ls_data'
 *
 */

/* Include files */
#include "mapTF2_ls_data.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
unsigned int state[625];

double freq;

bool freq_not_empty;

omp_nest_lock_t mapTF2_ls_nestLockGlobal;

const bool bv[128] = {
    false, false, false, false, false, false, false, false, false, true,  true,
    true,  true,  true,  false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, true,  true,  true,  true,  true,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false};

bool isInitialized_mapTF2_ls = false;

/* End of code generation (mapTF2_ls_data.c) */
