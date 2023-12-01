/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * schur.c
 *
 * Code generation for function 'schur'
 *
 */

/* Include files */
#include "schur.h"
#include "gkmPWM_rtwutil.h"
#include "xhseqr.h"
#include "xnrm2.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void schur(double A[16], double V[16])
{
  double work[4];
  double tau[3];
  double alpha1_tmp;
  double beta1;
  double temp;
  int alpha1_tmp_tmp;
  int b_i;
  int c_i;
  int coltop;
  int d_i;
  int exitg1;
  int i;
  int ia;
  int in;
  int jA;
  int knt;
  int lastc;
  int lastv;
  bool exitg2;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 0; i < 3; i++) {
    b_i = i << 2;
    in = (i + 1) << 2;
    alpha1_tmp_tmp = (i + b_i) + 1;
    alpha1_tmp = A[alpha1_tmp_tmp];
    if (i + 3 <= 4) {
      c_i = i + 1;
    } else {
      c_i = 2;
    }
    b_i = (c_i + b_i) + 2;
    tau[i] = 0.0;
    temp = xnrm2(2 - i, A, b_i);
    if (temp != 0.0) {
      beta1 = rt_hypotd(alpha1_tmp, temp);
      if (alpha1_tmp >= 0.0) {
        beta1 = -beta1;
      }
      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = 0;
        c_i = (b_i - i) + 1;
        do {
          knt++;
          for (coltop = b_i; coltop <= c_i; coltop++) {
            A[coltop - 1] *= 9.9792015476736E+291;
          }
          beta1 *= 9.9792015476736E+291;
          alpha1_tmp *= 9.9792015476736E+291;
        } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt < 20));
        beta1 = rt_hypotd(alpha1_tmp, xnrm2(2 - i, A, b_i));
        if (alpha1_tmp >= 0.0) {
          beta1 = -beta1;
        }
        tau[i] = (beta1 - alpha1_tmp) / beta1;
        temp = 1.0 / (alpha1_tmp - beta1);
        c_i = (b_i - i) + 1;
        for (coltop = b_i; coltop <= c_i; coltop++) {
          A[coltop - 1] *= temp;
        }
        for (coltop = 0; coltop < knt; coltop++) {
          beta1 *= 1.0020841800044864E-292;
        }
        alpha1_tmp = beta1;
      } else {
        tau[i] = (beta1 - alpha1_tmp) / beta1;
        temp = 1.0 / (alpha1_tmp - beta1);
        c_i = (b_i - i) + 1;
        for (coltop = b_i; coltop <= c_i; coltop++) {
          A[coltop - 1] *= temp;
        }
        alpha1_tmp = beta1;
      }
    }
    A[alpha1_tmp_tmp] = 1.0;
    jA = in + 1;
    if (tau[i] != 0.0) {
      lastv = 2 - i;
      b_i = (alpha1_tmp_tmp - i) + 2;
      while ((lastv + 1 > 0) && (A[b_i] == 0.0)) {
        lastv--;
        b_i--;
      }
      lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        knt = in + lastc;
        ia = knt;
        do {
          exitg1 = 0;
          if (ia <= knt + (lastv << 2)) {
            if (A[ia - 1] != 0.0) {
              exitg1 = 1;
            } else {
              ia += 4;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);
        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = -1;
      lastc = 0;
    }
    if (lastv + 1 > 0) {
      if (lastc != 0) {
        if (0 <= lastc - 1) {
          memset(&work[0], 0, lastc * sizeof(double));
        }
        knt = alpha1_tmp_tmp;
        c_i = (in + (lastv << 2)) + 1;
        for (b_i = jA; b_i <= c_i; b_i += 4) {
          d_i = (b_i + lastc) - 1;
          for (ia = b_i; ia <= d_i; ia++) {
            coltop = ia - b_i;
            work[coltop] += A[ia - 1] * A[knt];
          }
          knt++;
        }
      }
      if (-tau[i] != 0.0) {
        jA = in;
        for (knt = 0; knt <= lastv; knt++) {
          temp = A[alpha1_tmp_tmp + knt];
          if (temp != 0.0) {
            temp *= -tau[i];
            c_i = jA + 1;
            d_i = lastc + jA;
            for (b_i = c_i; b_i <= d_i; b_i++) {
              A[b_i - 1] += work[(b_i - jA) - 1] * temp;
            }
          }
          jA += 4;
        }
      }
    }
    jA = (i + in) + 2;
    if (tau[i] != 0.0) {
      lastv = 3 - i;
      b_i = (alpha1_tmp_tmp - i) + 2;
      while ((lastv > 0) && (A[b_i] == 0.0)) {
        lastv--;
        b_i--;
      }
      lastc = 2 - i;
      exitg2 = false;
      while ((!exitg2) && (lastc + 1 > 0)) {
        coltop = jA + (lastc << 2);
        ia = coltop;
        do {
          exitg1 = 0;
          if (ia <= (coltop + lastv) - 1) {
            if (A[ia - 1] != 0.0) {
              exitg1 = 1;
            } else {
              ia++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);
        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = -1;
    }
    if (lastv > 0) {
      if (lastc + 1 != 0) {
        if (0 <= lastc) {
          memset(&work[0], 0, (lastc + 1) * sizeof(double));
        }
        c_i = jA + (lastc << 2);
        for (b_i = jA; b_i <= c_i; b_i += 4) {
          temp = 0.0;
          d_i = (b_i + lastv) - 1;
          for (ia = b_i; ia <= d_i; ia++) {
            temp += A[ia - 1] * A[(alpha1_tmp_tmp + ia) - b_i];
          }
          coltop = (b_i - jA) >> 2;
          work[coltop] += temp;
        }
      }
      if (-tau[i] != 0.0) {
        for (knt = 0; knt <= lastc; knt++) {
          temp = work[knt];
          if (temp != 0.0) {
            temp *= -tau[i];
            c_i = lastv + jA;
            for (b_i = jA; b_i < c_i; b_i++) {
              A[b_i - 1] += A[(alpha1_tmp_tmp + b_i) - jA] * temp;
            }
          }
          jA += 4;
        }
      }
    }
    A[alpha1_tmp_tmp] = alpha1_tmp;
  }
  memcpy(&V[0], &A[0], 16U * sizeof(double));
  xhseqr(V);
}

/* End of code generation (schur.c) */
