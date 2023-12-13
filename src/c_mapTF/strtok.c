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
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_strtok(const emxArray_char_T *x, emxArray_char_T *token,
              emxArray_char_T *remain)
{
  static const char b_cv[6] = {'\x09', '\x0a', '\x0b', '\x0c', '\x0d', ' '};
  int b_k;
  int exitg1;
  int i;
  int itoken;
  int k;
  int loop_ub;
  int n;
  const char *x_data;
  char *remain_data;
  bool exitg2;
  x_data = x->data;
  n = x->size[1];
  k = 0;
  exitg2 = false;
  while ((!exitg2) && (k + 1 <= n)) {
    b_k = 0;
    do {
      exitg1 = 0;
      if (b_k < 6) {
        if (x_data[k] == b_cv[b_k]) {
          k++;
          exitg1 = 1;
        } else {
          b_k++;
        }
      } else {
        exitg1 = 2;
      }
    } while (exitg1 == 0);
    if (exitg1 != 1) {
      exitg2 = true;
    }
  }
  itoken = k + 1;
  exitg2 = false;
  while ((!exitg2) && (k + 1 <= n)) {
    b_k = 0;
    do {
      exitg1 = 0;
      if (b_k < 6) {
        if (x_data[k] == b_cv[b_k]) {
          exitg1 = 1;
        } else {
          b_k++;
        }
      } else {
        k++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);
    if (exitg1 == 1) {
      exitg2 = true;
    }
  }
  if (k + 1 > x->size[1]) {
    n = 0;
    i = 0;
  } else {
    n = k;
    i = x->size[1];
  }
  b_k = remain->size[0] * remain->size[1];
  remain->size[0] = 1;
  loop_ub = i - n;
  remain->size[1] = loop_ub;
  emxEnsureCapacity_char_T(remain, b_k);
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
void d_strtok(const emxArray_char_T *x, emxArray_char_T *token,
              emxArray_char_T *remain)
{
  static const char b_cv[2] = {' ', '\x09'};
  int b_k;
  int exitg1;
  int i;
  int itoken;
  int k;
  int loop_ub;
  int n;
  const char *x_data;
  char *remain_data;
  bool exitg2;
  x_data = x->data;
  n = x->size[1];
  k = 0;
  exitg2 = false;
  while ((!exitg2) && (k + 1 <= n)) {
    b_k = 0;
    do {
      exitg1 = 0;
      if (b_k < 2) {
        if (x_data[k] == b_cv[b_k]) {
          k++;
          exitg1 = 1;
        } else {
          b_k++;
        }
      } else {
        exitg1 = 2;
      }
    } while (exitg1 == 0);
    if (exitg1 != 1) {
      exitg2 = true;
    }
  }
  itoken = k + 1;
  exitg2 = false;
  while ((!exitg2) && (k + 1 <= n)) {
    b_k = 0;
    do {
      exitg1 = 0;
      if (b_k < 2) {
        if (x_data[k] == b_cv[b_k]) {
          exitg1 = 1;
        } else {
          b_k++;
        }
      } else {
        k++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);
    if (exitg1 == 1) {
      exitg2 = true;
    }
  }
  if (k + 1 > x->size[1]) {
    n = 0;
    i = 0;
  } else {
    n = k;
    i = x->size[1];
  }
  b_k = remain->size[0] * remain->size[1];
  remain->size[0] = 1;
  loop_ub = i - n;
  remain->size[1] = loop_ub;
  emxEnsureCapacity_char_T(remain, b_k);
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
