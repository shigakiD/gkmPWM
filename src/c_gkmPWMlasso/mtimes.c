/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.c
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include "cblas.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_mtimes(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *C)
{
  const double *A_data;
  const double *B_data;
  double *C_data;
  int i;
  int loop_ub;
  B_data = B->data;
  A_data = A->data;
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0) ||
      (B->size[1] == 0)) {
    i = C->size[0] * C->size[1];
    C->size[0] = A->size[1];
    C->size[1] = B->size[1];
    emxEnsureCapacity_real_T(C, i);
    C_data = C->data;
    loop_ub = A->size[1] * B->size[1];
    for (i = 0; i < loop_ub; i++) {
      C_data[i] = 0.0;
    }
  } else {
    i = C->size[0] * C->size[1];
    C->size[0] = A->size[1];
    C->size[1] = B->size[1];
    emxEnsureCapacity_real_T(C, i);
    C_data = C->data;
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)A->size[1],
                (blasint)B->size[1], (blasint)A->size[0], 1.0, &A_data[0],
                (blasint)A->size[0], &B_data[0], (blasint)B->size[0], 0.0,
                &C_data[0], (blasint)A->size[1]);
  }
}

/*
 *
 */
void mtimes(const emxArray_real_T *A, const emxArray_real_T *B,
            emxArray_real_T *C)
{
  const double *A_data;
  const double *B_data;
  double *C_data;
  int i;
  int loop_ub;
  B_data = B->data;
  A_data = A->data;
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0) ||
      (B->size[1] == 0)) {
    i = C->size[0] * C->size[1];
    C->size[0] = A->size[0];
    C->size[1] = B->size[0];
    emxEnsureCapacity_real_T(C, i);
    C_data = C->data;
    loop_ub = A->size[0] * B->size[0];
    for (i = 0; i < loop_ub; i++) {
      C_data[i] = 0.0;
    }
  } else {
    i = C->size[0] * C->size[1];
    C->size[0] = A->size[0];
    C->size[1] = B->size[0];
    emxEnsureCapacity_real_T(C, i);
    C_data = C->data;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)A->size[0],
                (blasint)B->size[0], (blasint)A->size[1], 1.0, &A_data[0],
                (blasint)A->size[0], &B_data[0], (blasint)B->size[0], 0.0,
                &C_data[0], (blasint)A->size[0]);
  }
}

/* End of code generation (mtimes.c) */
