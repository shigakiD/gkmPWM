/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * print_processing.c
 *
 * Code generation for function 'print_processing'
 *
 */

/* Include files */
#include "print_processing.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void print_processing(const emxArray_char_T *varargin_1,
                      const emxArray_char_T *varargin_2,
                      const emxArray_char_T *varargin_3,
                      cell_wrap_8 validatedArgumetns[3])
{
  int i;
  int loop_ub;
  const char *varargin_1_data;
  const char *varargin_2_data;
  const char *varargin_3_data;
  varargin_3_data = varargin_3->data;
  varargin_2_data = varargin_2->data;
  varargin_1_data = varargin_1->data;
  i = validatedArgumetns[0].f1->size[0] * validatedArgumetns[0].f1->size[1];
  validatedArgumetns[0].f1->size[0] = 1;
  validatedArgumetns[0].f1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(validatedArgumetns[0].f1, i);
  loop_ub = varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    validatedArgumetns[0].f1->data[i] = varargin_1_data[i];
  }
  validatedArgumetns[0].f1->data[varargin_1->size[1]] = '\x00';
  i = validatedArgumetns[1].f1->size[0] * validatedArgumetns[1].f1->size[1];
  validatedArgumetns[1].f1->size[0] = 1;
  validatedArgumetns[1].f1->size[1] = varargin_2->size[1] + 1;
  emxEnsureCapacity_char_T(validatedArgumetns[1].f1, i);
  loop_ub = varargin_2->size[1];
  for (i = 0; i < loop_ub; i++) {
    validatedArgumetns[1].f1->data[i] = varargin_2_data[i];
  }
  validatedArgumetns[1].f1->data[varargin_2->size[1]] = '\x00';
  i = validatedArgumetns[2].f1->size[0] * validatedArgumetns[2].f1->size[1];
  validatedArgumetns[2].f1->size[0] = 1;
  validatedArgumetns[2].f1->size[1] = varargin_3->size[1] + 1;
  emxEnsureCapacity_char_T(validatedArgumetns[2].f1, i);
  loop_ub = varargin_3->size[1];
  for (i = 0; i < loop_ub; i++) {
    validatedArgumetns[2].f1->data[i] = varargin_3_data[i];
  }
  validatedArgumetns[2].f1->data[varargin_3->size[1]] = '\x00';
}

/* End of code generation (print_processing.c) */
