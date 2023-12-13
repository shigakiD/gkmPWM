/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * isequal.c
 *
 * Code generation for function 'isequal'
 *
 */

/* Include files */
#include "isequal.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
bool isequal(const emxArray_boolean_T *varargin_1,
             const emxArray_boolean_T *varargin_2)
{
  int k;
  const bool *varargin_1_data;
  const bool *varargin_2_data;
  bool exitg1;
  bool p;
  varargin_2_data = varargin_2->data;
  varargin_1_data = varargin_1->data;
  p = (varargin_1->size[1] == varargin_2->size[1]);
  if (p && (varargin_1->size[1] != 0) && (varargin_2->size[1] != 0)) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= varargin_2->size[1] - 1)) {
      if (varargin_1_data[k] != varargin_2_data[k]) {
        p = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }
  return p;
}

/* End of code generation (isequal.c) */
