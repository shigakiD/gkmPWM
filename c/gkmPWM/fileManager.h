/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fileManager.h
 *
 * Code generation for function 'fileManager'
 *
 */

#ifndef FILEMANAGER_H
#define FILEMANAGER_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
int cfclose(double fid);

signed char cfopen(const emxArray_char_T *cfilename, const char *cpermission);

void filedata_init(void);

void getfilestar(double fid, FILE **filestar, bool *autoflush);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (fileManager.h) */
