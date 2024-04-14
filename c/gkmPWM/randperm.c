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
  c_rand(p);
  for (k = 0; k <= 2; k += 2) {
    if (p[k] <= p[k + 1]) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
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
  for (i = 0; i < 4; i++) {
    p[i] = idx[i];
  }
}

/* End of code generation (randperm.c) */
