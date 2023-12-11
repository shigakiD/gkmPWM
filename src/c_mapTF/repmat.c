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
#include "mapTF2_ls_emxutil.h"
#include "mapTF2_ls_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void repmat(const double a[4], double varargin_1, emxArray_real_T *b)
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

/* End of code generation (repmat.c) */
