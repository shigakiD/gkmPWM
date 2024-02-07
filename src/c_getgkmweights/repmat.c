/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "repmat.h"

/* Function Definitions */
/*
 *
 */
void repmat(const double a[4], double b[16])
{
  double d;
  double d1;
  double d2;
  double d3;
  int ibtile;
  int jtilecol;
  d = a[0];
  d1 = a[1];
  d2 = a[2];
  d3 = a[3];
  for (jtilecol = 0; jtilecol < 4; jtilecol++) {
    ibtile = jtilecol << 2;
    b[ibtile] = d;
    b[ibtile + 1] = d1;
    b[ibtile + 2] = d2;
    b[ibtile + 3] = d3;
  }
}

/* End of code generation (repmat.c) */
