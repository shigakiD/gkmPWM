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
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_sprintf(const emxArray_char_T *varargin_1, int varargin_2,
               int varargin_3, emxArray_char_T *str_Value)
{
  emxArray_char_T *b_varargin_1;
  emxArray_char_T *c_varargin_1;
  int i;
  int nbytes;
  const char *varargin_1_data;
  char *b_varargin_1_data;
  char *c_varargin_1_data;
  varargin_1_data = varargin_1->data;
  emxInit_char_T(&b_varargin_1, 2);
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  nbytes = varargin_1->size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_1_data[i] = varargin_1_data[i];
  }
  emxInit_char_T(&c_varargin_1, 2);
  b_varargin_1_data[varargin_1->size[1]] = '\x00';
  i = c_varargin_1->size[0] * c_varargin_1->size[1];
  c_varargin_1->size[0] = 1;
  c_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_1, i);
  c_varargin_1_data = c_varargin_1->data;
  nbytes = varargin_1->size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_1_data[i] = varargin_1_data[i];
  }
  c_varargin_1_data[varargin_1->size[1]] = '\x00';
  nbytes = snprintf(NULL, 0, "%s_%d_%d", &c_varargin_1_data[0], varargin_2,
                    varargin_3);
  i = str_Value->size[0] * str_Value->size[1];
  str_Value->size[0] = 1;
  str_Value->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(str_Value, i);
  c_varargin_1_data = str_Value->data;
  snprintf(&c_varargin_1_data[0], (size_t)(nbytes + 1), "%s_%d_%d",
           &b_varargin_1_data[0], varargin_2, varargin_3);
  i = str_Value->size[0] * str_Value->size[1];
  if (1 > nbytes) {
    str_Value->size[1] = 0;
  } else {
    str_Value->size[1] = nbytes;
  }
  emxEnsureCapacity_char_T(str_Value, i);
  emxFree_char_T(&b_varargin_1);
  emxFree_char_T(&c_varargin_1);
}

/* End of code generation (sprintf.c) */
