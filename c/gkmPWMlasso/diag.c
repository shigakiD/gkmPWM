/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diag.c
 *
 * Code generation for function 'diag'
 *
 */

/* Include files */
#include "diag.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void diag(const emxArray_real_T *v, emxArray_real_T *d)
{
  const double *v_data;
  double *d_data;
  int i;
  int loop_ub;
  int nv;
  v_data = v->data;
  nv = v->size[0];
  i = d->size[0] * d->size[1];
  d->size[0] = v->size[0];
  d->size[1] = v->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  loop_ub = v->size[0] * v->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = 0.0;
  }
  for (loop_ub = 0; loop_ub < nv; loop_ub++) {
    d_data[loop_ub + d->size[0] * loop_ub] = v_data[loop_ub];
  }
}

/* End of code generation (diag.c) */
