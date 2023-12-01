/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sprintf.c
 *
 * Code generation for function 'sprintf'
 *
 */

/* Include files */
#include "sprintf.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_sprintf(const emxArray_char_T *varargin_1_Value, emxArray_char_T *str)
{
  emxArray_char_T *b_varargin_1;
  emxArray_char_T *varargin_1;
  int i;
  int nbytes;
  const char *varargin_1_Value_data;
  char *str_data;
  char *varargin_1_data;
  varargin_1_Value_data = varargin_1_Value->data;
  emxInit_char_T(&varargin_1, 2);
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = 1;
  varargin_1->size[1] = varargin_1_Value->size[1] + 1;
  emxEnsureCapacity_char_T(varargin_1, i);
  varargin_1_data = varargin_1->data;
  nbytes = varargin_1_Value->size[1];
  for (i = 0; i < nbytes; i++) {
    varargin_1_data[i] = varargin_1_Value_data[i];
  }
  emxInit_char_T(&b_varargin_1, 2);
  varargin_1_data[varargin_1_Value->size[1]] = '\x00';
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1_Value->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  str_data = b_varargin_1->data;
  nbytes = varargin_1_Value->size[1];
  for (i = 0; i < nbytes; i++) {
    str_data[i] = varargin_1_Value_data[i];
  }
  str_data[varargin_1_Value->size[1]] = '\x00';
  nbytes = snprintf(NULL, 0, "%s_gkmPWM.out", &str_data[0]);
  i = str->size[0] * str->size[1];
  str->size[0] = 1;
  str->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(str, i);
  str_data = str->data;
  snprintf(&str_data[0], (size_t)(nbytes + 1), "%s_gkmPWM.out",
           &varargin_1_data[0]);
  i = str->size[0] * str->size[1];
  if (1 > nbytes) {
    str->size[1] = 0;
  } else {
    str->size[1] = nbytes;
  }
  emxEnsureCapacity_char_T(str, i);
  emxFree_char_T(&varargin_1);
  emxFree_char_T(&b_varargin_1);
}

/* End of code generation (sprintf.c) */
