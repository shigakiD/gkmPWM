/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * strip.c
 *
 * Code generation for function 'strip'
 *
 */

/* Include files */
#include "strip.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void strip(const emxArray_char_T *str, emxArray_char_T *s)
{
  static const char stripCharacters[10] = {'\x09', '\x0a', '\x0b', '\x0c',
                                           '\x0d', '\x1c', '\x1d', '\x1e',
                                           '\x1f', ' '};
  int b_i;
  int beginningIndex;
  int endIndex;
  int i;
  int loop_ub;
  const char *str_data;
  char b_str;
  char *s_data;
  bool x[10];
  bool exitg1;
  bool guard1 = false;
  bool nonStripCharFound;
  bool y;
  str_data = str->data;
  if (str->size[1] == 0) {
    s->size[0] = 1;
    s->size[1] = 0;
  } else {
    i = 1;
    beginningIndex = 0;
    nonStripCharFound = false;
    while (!nonStripCharFound) {
      guard1 = false;
      if (i == str->size[1] + 1) {
        guard1 = true;
      } else {
        b_str = str_data[i - 1];
        for (b_i = 0; b_i < 10; b_i++) {
          x[b_i] = (b_str != stripCharacters[b_i]);
        }
        y = true;
        loop_ub = 0;
        exitg1 = false;
        while ((!exitg1) && (loop_ub < 10)) {
          if (!x[loop_ub]) {
            y = false;
            exitg1 = true;
          } else {
            loop_ub++;
          }
        }
        if (y) {
          guard1 = true;
        }
      }
      if (guard1) {
        nonStripCharFound = true;
        beginningIndex = i;
      }
      i++;
    }
    i = str->size[1];
    endIndex = 0;
    nonStripCharFound = false;
    while (!nonStripCharFound) {
      guard1 = false;
      if (i == 0) {
        guard1 = true;
      } else {
        b_str = str_data[i - 1];
        for (b_i = 0; b_i < 10; b_i++) {
          x[b_i] = (b_str != stripCharacters[b_i]);
        }
        y = true;
        loop_ub = 0;
        exitg1 = false;
        while ((!exitg1) && (loop_ub < 10)) {
          if (!x[loop_ub]) {
            y = false;
            exitg1 = true;
          } else {
            loop_ub++;
          }
        }
        if (y) {
          guard1 = true;
        }
      }
      if (guard1) {
        nonStripCharFound = true;
        endIndex = i;
      }
      i--;
    }
    if (beginningIndex > endIndex) {
      b_i = 0;
      endIndex = 0;
    } else {
      b_i = beginningIndex - 1;
    }
    i = s->size[0] * s->size[1];
    s->size[0] = 1;
    loop_ub = endIndex - b_i;
    s->size[1] = loop_ub;
    emxEnsureCapacity_char_T(s, i);
    s_data = s->data;
    for (i = 0; i < loop_ub; i++) {
      s_data[i] = str_data[b_i + i];
    }
  }
}

/* End of code generation (strip.c) */
