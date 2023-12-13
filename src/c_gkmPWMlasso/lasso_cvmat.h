/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * lasso_cvmat.h
 *
 * Code generation for function 'lasso_cvmat'
 *
 */

#ifndef LASSO_CVMAT_H
#define LASSO_CVMAT_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_lasso_cvmat(emxArray_real_T *X, emxArray_real_T *Y, double dfMax,
                   emxArray_real_T *B, double stats_Intercept_data[],
                   int stats_Intercept_size[2], double stats_Lambda_data[],
                   int stats_Lambda_size[2], double *stats_Alpha,
                   double stats_DF_data[], int stats_DF_size[2],
                   double stats_MSE_data[], int stats_MSE_size[2]);

void lasso_cvmat(emxArray_real_T *X, emxArray_real_T *Y, double dfMax,
                 emxArray_real_T *B, double stats_Intercept_data[],
                 int stats_Intercept_size[2], double stats_Lambda_data[],
                 int stats_Lambda_size[2], double *stats_Alpha,
                 double stats_DF_data[], int stats_DF_size[2],
                 double stats_MSE_data[], int stats_MSE_size[2]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (lasso_cvmat.h) */
