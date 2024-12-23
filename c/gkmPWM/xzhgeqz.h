/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.h
 *
 * Code generation for function 'xzhgeqz'
 *
 */

#ifndef XZHGEQZ_H
#define XZHGEQZ_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void xzhgeqz(const creal_T A[16], int ilo, int ihi, int *info,
             creal_T alpha1[4], creal_T beta1[4]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xzhgeqz.h) */
