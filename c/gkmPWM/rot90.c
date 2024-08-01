/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rot90.c
 *
 * Code generation for function 'rot90'
 *
 */

/* Include files */
#include "rot90.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_rot90(const double A[16], double B[16])
{
  int i;
  int j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      B[i + (j << 2)] = A[((i << 2) - j) + 3];
    }
  }
}

/*
 *
 */
void c_rot90(const double A[16], double B[16])
{
  int i;
  int j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      B[i + (j << 2)] = A[j + ((3 - i) << 2)];
    }
  }
}

/*
 *
 */
void d_rot90(const emxArray_real_T *A, emxArray_real_T *B)
{
  const double *A_data;
  double *B_data;
  int b_i;
  int b_j;
  int i;
  int j;
  int m;
  A_data = A->data;
  m = A->size[0];
  j = B->size[0] * B->size[1];
  B->size[0] = A->size[0];
  B->size[1] = 4;
  emxEnsureCapacity_real_T(B, j);
  B_data = B->data;
  if ((A->size[0] << 2) >= 8192) {

    for (b_j = 0; b_j < 4; b_j++) {
      for (b_i = 0; b_i < m; b_i++) {
        B_data[b_i + B->size[0] * b_j] =
            A_data[((m - b_i) + A->size[0] * (3 - b_j)) - 1];
      }
    }
  } else {
    for (j = 0; j < 4; j++) {
      for (i = 0; i < m; i++) {
        B_data[i + B->size[0] * j] =
            A_data[((m - i) + A->size[0] * (3 - j)) - 1];
      }
    }
  }
}

/*
 *
 */
void rot90(const emxArray_real_T *A, emxArray_real_T *B)
{
  const double *A_data;
  double *B_data;
  int b_i;
  int b_j;
  int i;
  int j;
  int m;
  int n;
  A_data = A->data;
  m = A->size[0];
  n = A->size[1];
  j = B->size[0] * B->size[1];
  B->size[0] = A->size[0];
  B->size[1] = A->size[1];
  emxEnsureCapacity_real_T(B, j);
  B_data = B->data;
  if (A->size[0] * A->size[1] >= 8192) {
    j = A->size[1] - 1;

    for (b_j = 0; b_j <= j; b_j++) {
      for (b_i = 0; b_i < m; b_i++) {
        B_data[b_i + B->size[0] * b_j] =
            A_data[((m - b_i) + A->size[0] * ((n - b_j) - 1)) - 1];
      }
    }
  } else {
    for (j = 0; j < n; j++) {
      for (i = 0; i < m; i++) {
        B_data[i + B->size[0] * j] =
            A_data[((m - i) + A->size[0] * ((n - j) - 1)) - 1];
      }
    }
  }
}

/* End of code generation (rot90.c) */
