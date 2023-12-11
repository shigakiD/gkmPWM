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
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_repmat(const double a[4], double varargin_1, emxArray_real_T *b)
{
  double *b_data;
  int i;
  int ibmat;
  int itilerow;
  int jcol;
  i = (int)varargin_1;
  ibmat = b->size[0] * b->size[1];
  b->size[0] = (int)varargin_1;
  b->size[1] = 4;
  emxEnsureCapacity_real_T(b, ibmat);
  b_data = b->data;
  for (jcol = 0; jcol < 4; jcol++) {
    ibmat = jcol * (int)varargin_1;
    for (itilerow = 0; itilerow < i; itilerow++) {
      b_data[ibmat + itilerow] = a[jcol];
    }
  }
}

/*
 *
 */
void c_repmat(const emxArray_real_T *a, emxArray_real_T *b)
{
  const double *a_data;
  double *b_data;
  int ibcol;
  int itilerow;
  int k;
  int nrows;
  a_data = a->data;
  nrows = b->size[0];
  b->size[0] = a->size[0] << 2;
  emxEnsureCapacity_real_T(b, nrows);
  b_data = b->data;
  nrows = a->size[0];
  for (itilerow = 0; itilerow < 4; itilerow++) {
    ibcol = itilerow * nrows;
    for (k = 0; k < nrows; k++) {
      b_data[ibcol + k] = a_data[k];
    }
  }
}

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
