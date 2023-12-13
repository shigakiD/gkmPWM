/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fgetl.c
 *
 * Code generation for function 'fgetl'
 *
 */

/* Include files */
#include "fgetl.h"
#include "fgets.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void fgetl(double fileID, emxArray_char_T *out)
{
  int i;
  int loop_ub;
  char *out_data;
  b_fgets(fileID, out);
  out_data = out->data;
  if (out->size[1] != 0) {
    if (out_data[out->size[1] - 1] == '\x0a') {
      if ((out->size[1] > 1) && (out_data[out->size[1] - 2] == '\x0d')) {
        if (1 > out->size[1] - 2) {
          loop_ub = 0;
        } else {
          loop_ub = out->size[1] - 2;
        }
        i = out->size[0] * out->size[1];
        out->size[0] = 1;
        out->size[1] = loop_ub;
        emxEnsureCapacity_char_T(out, i);
      } else {
        if (1 > out->size[1] - 1) {
          loop_ub = 0;
        } else {
          loop_ub = out->size[1] - 1;
        }
        i = out->size[0] * out->size[1];
        out->size[0] = 1;
        out->size[1] = loop_ub;
        emxEnsureCapacity_char_T(out, i);
      }
    } else if (out_data[out->size[1] - 1] == '\x0d') {
      if (1 > out->size[1] - 1) {
        loop_ub = 0;
      } else {
        loop_ub = out->size[1] - 1;
      }
      i = out->size[0] * out->size[1];
      out->size[0] = 1;
      out->size[1] = loop_ub;
      emxEnsureCapacity_char_T(out, i);
    }
  }
}

/* End of code generation (fgetl.c) */
