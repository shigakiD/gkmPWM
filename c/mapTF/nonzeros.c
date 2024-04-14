/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nonzeros.c
 *
 * Code generation for function 'nonzeros'
 *
 */

/* Include files */
#include "nonzeros.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void nonzeros(const emxArray_real_T *s, emxArray_real_T *v)
{
  const double *s_data;
  double *v_data;
  int b_i;
  int i;
  int k;
  int n;
  s_data = s->data;
  n = (s->size[1] << 2) * s->size[2];
  i = 0;
  b_i = (s->size[1] << 2) * s->size[2];
  for (k = 0; k < b_i; k++) {
    if (s_data[k] != 0.0) {
      i++;
    }
  }
  b_i = v->size[0];
  v->size[0] = i;
  emxEnsureCapacity_real_T(v, b_i);
  v_data = v->data;
  i = -1;
  for (k = 0; k < n; k++) {
    if (s_data[k] != 0.0) {
      i++;
      v_data[i] = s_data[k];
    }
  }
}

/* End of code generation (nonzeros.c) */
