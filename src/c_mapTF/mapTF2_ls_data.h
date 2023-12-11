/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF2_ls_data.h
 *
 * Code generation for function 'mapTF2_ls_data'
 *
 */

#ifndef MAPTF2_LS_DATA_H
#define MAPTF2_LS_DATA_H

/* Include files */
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

/* Variable Declarations */
extern unsigned int state[625];
extern double freq;
extern bool freq_not_empty;
extern omp_nest_lock_t mapTF2_ls_nestLockGlobal;
extern const bool bv[128];
extern bool isInitialized_mapTF2_ls;

#endif
/* End of code generation (mapTF2_ls_data.h) */
