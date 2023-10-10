/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * strtok.c
 *
 * Code generation for function 'strtok'
 *
 */

/* Include files */
#include "strtok.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_strtok(const emxArray_char_T *x, emxArray_char_T *token,
              emxArray_char_T *remain)
{
  int i;
  int i1;
  int itoken;
  int k;
  int loop_ub;
  int n;
  const char *x_data;
  char *remain_data;
  x_data = x->data;
  n = x->size[1];
  k = 0;
  while ((k + 1 <= n) && (x_data[k] == '\x09')) {
    k++;
  }
  itoken = k + 1;
  while ((k + 1 <= n) && (!(x_data[k] == '\x09'))) {
    k++;
  }
  if (k + 1 > x->size[1]) {
    n = 0;
    i = 0;
  } else {
    n = k;
    i = x->size[1];
  }
  i1 = remain->size[0] * remain->size[1];
  remain->size[0] = 1;
  loop_ub = i - n;
  remain->size[1] = loop_ub;
  emxEnsureCapacity_char_T(remain, i1);
  remain_data = remain->data;
  for (i = 0; i < loop_ub; i++) {
    remain_data[i] = x_data[n + i];
  }
  if (itoken > k) {
    n = 0;
    k = 0;
  } else {
    n = itoken - 1;
  }
  i = token->size[0] * token->size[1];
  token->size[0] = 1;
  loop_ub = k - n;
  token->size[1] = loop_ub;
  emxEnsureCapacity_char_T(token, i);
  remain_data = token->data;
  for (i = 0; i < loop_ub; i++) {
    remain_data[i] = x_data[n + i];
  }
}

/*
 *
 */
void c_strtok(const emxArray_char_T *x, emxArray_char_T *token,
              emxArray_char_T *remain)
{
  int i;
  int i1;
  int itoken;
  int k;
  int loop_ub;
  int n;
  const char *x_data;
  char *remain_data;
  x_data = x->data;
  n = x->size[1];
  k = 0;
  while ((k + 1 <= n) && (x_data[k] == ' ')) {
    k++;
  }
  itoken = k + 1;
  while ((k + 1 <= n) && (!(x_data[k] == ' '))) {
    k++;
  }
  if (k + 1 > x->size[1]) {
    n = 0;
    i = 0;
  } else {
    n = k;
    i = x->size[1];
  }
  i1 = remain->size[0] * remain->size[1];
  remain->size[0] = 1;
  loop_ub = i - n;
  remain->size[1] = loop_ub;
  emxEnsureCapacity_char_T(remain, i1);
  remain_data = remain->data;
  for (i = 0; i < loop_ub; i++) {
    remain_data[i] = x_data[n + i];
  }
  if (itoken > k) {
    n = 0;
    k = 0;
  } else {
    n = itoken - 1;
  }
  i = token->size[0] * token->size[1];
  token->size[0] = 1;
  loop_ub = k - n;
  token->size[1] = loop_ub;
  emxEnsureCapacity_char_T(token, i);
  remain_data = token->data;
  for (i = 0; i < loop_ub; i++) {
    remain_data[i] = x_data[n + i];
  }
}

/* End of code generation (strtok.c) */
