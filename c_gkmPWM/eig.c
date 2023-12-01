/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eig.c
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include "eig.h"
#include "schur.h"
#include "lapacke.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void eig(const double A[16], creal_T V[4])
{
  lapack_int ihi_t;
  lapack_int ilo_t;
  lapack_int info_t;
  double T[16];
  double b_A[16];
  double scale[4];
  double wimag[4];
  double wreal[4];
  double abnrm;
  double rconde;
  double rcondv;
  double vleft;
  double vright;
  int exitg1;
  int i;
  int j;
  bool exitg2;
  bool guard1 = false;
  bool p;
  p = true;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 4)) {
    i = 0;
    do {
      exitg1 = 0;
      if (i <= j) {
        if (A[i + (j << 2)] != A[j + (i << 2)]) {
          p = false;
          exitg1 = 1;
        } else {
          i++;
        }
      } else {
        j++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);
    if (exitg1 == 1) {
      exitg2 = true;
    }
  }
  if (p) {
    memcpy(&b_A[0], &A[0], 16U * sizeof(double));
    schur(b_A, T);
    V[0].re = T[0];
    V[0].im = 0.0;
    V[1].re = T[5];
    V[1].im = 0.0;
    V[2].re = T[10];
    V[2].im = 0.0;
    V[3].re = T[15];
    V[3].im = 0.0;
  } else {
    p = true;
    j = 0;
    exitg2 = false;
    while ((!exitg2) && (j < 4)) {
      i = 0;
      do {
        exitg1 = 0;
        if (i <= j) {
          if (A[i + (j << 2)] != -A[j + (i << 2)]) {
            p = false;
            exitg1 = 1;
          } else {
            i++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
    if (p) {
      memcpy(&b_A[0], &A[0], 16U * sizeof(double));
      schur(b_A, T);
      i = 1;
      do {
        exitg1 = 0;
        if (i <= 4) {
          guard1 = false;
          if (i != 4) {
            vleft = T[i + ((i - 1) << 2)];
            if (vleft != 0.0) {
              vleft = fabs(vleft);
              V[i - 1].re = 0.0;
              V[i - 1].im = vleft;
              V[i].re = 0.0;
              V[i].im = -vleft;
              i += 2;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
          if (guard1) {
            V[i - 1].re = 0.0;
            V[i - 1].im = 0.0;
            i++;
          }
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    } else {
      memcpy(&T[0], &A[0], 16U * sizeof(double));
      info_t = LAPACKE_dgeevx(
          LAPACK_COL_MAJOR, 'B', 'N', 'N', 'N', (lapack_int)4, &T[0],
          (lapack_int)4, &wreal[0], &wimag[0], &vleft, (lapack_int)1, &vright,
          (lapack_int)1, &ilo_t, &ihi_t, &scale[0], &abnrm, &rconde, &rcondv);
      if ((int)info_t < 0) {
        memset(&V[0], 0, 4U * sizeof(creal_T));
      } else {
        V[0].re = wreal[0];
        V[0].im = wimag[0];
        V[1].re = wreal[1];
        V[1].im = wimag[1];
        V[2].re = wreal[2];
        V[2].im = wimag[2];
        V[3].re = wreal[3];
        V[3].im = wimag[3];
      }
    }
  }
}

/* End of code generation (eig.c) */
