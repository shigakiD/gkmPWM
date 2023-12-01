/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * randperm.c
 *
 * Code generation for function 'randperm'
 *
 */

/* Include files */
#include "randperm.h"
#include "rand.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void randperm(double p[4])
{
  int idx[4];
  int iwork[4];
  int b_p;
  int i;
  int i1;
  int j;
  int k;
  int pEnd;
  int q;
  int qEnd;
  b_rand(p);
  if (p[0] <= p[1]) {
    idx[0] = 1;
    idx[1] = 2;
  } else {
    idx[0] = 2;
    idx[1] = 1;
  }
  if (p[2] <= p[3]) {
    idx[2] = 3;
    idx[3] = 4;
  } else {
    idx[2] = 4;
    idx[3] = 3;
  }
  j = 1;
  for (pEnd = 3; pEnd < 5; pEnd = qEnd + 2) {
    b_p = j;
    q = pEnd;
    qEnd = j + 4;
    k = 0;
    while (k + 1 <= 4) {
      i = idx[q - 1];
      i1 = idx[b_p - 1];
      if (p[i1 - 1] <= p[i - 1]) {
        iwork[k] = i1;
        b_p++;
        if (b_p == pEnd) {
          while (q < j + 4) {
            k++;
            iwork[k] = idx[q - 1];
            q++;
          }
        }
      } else {
        iwork[k] = i;
        q++;
        if (q == j + 4) {
          while (b_p < pEnd) {
            k++;
            iwork[k] = idx[b_p - 1];
            b_p++;
          }
        }
      }
      k++;
    }
    for (k = 0; k < 4; k++) {
      idx[(j + k) - 1] = iwork[k];
    }
    j += 4;
  }
  p[0] = idx[0];
  p[1] = idx[1];
  p[2] = idx[2];
  p[3] = idx[3];
}

/* End of code generation (randperm.c) */
