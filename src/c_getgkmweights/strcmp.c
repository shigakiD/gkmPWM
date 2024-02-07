/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * strcmp.c
 *
 * Code generation for function 'strcmp'
 *
 */

/* Include files */
#include "strcmp.h"
#include "getgkmweights_types.h"

/* Function Definitions */
/*
 *
 */
bool b_strcmp(const emxArray_char_T *a)
{
  static const char cv[3] = {'a', 'l', 'l'};
  int exitg1;
  int kstr;
  const char *a_data;
  bool b_bool;
  a_data = a->data;
  b_bool = false;
  if (a->size[1] == 3) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 3) {
        if (a_data[kstr] != cv[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return b_bool;
}

/*
 *
 */
bool c_strcmp(const emxArray_char_T *b)
{
  static const char cv[2] = {'S', 'V'};
  int exitg1;
  int kstr;
  const char *b_data;
  bool b_bool;
  b_data = b->data;
  b_bool = false;
  if (2 == b->size[1]) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 2) {
        if (cv[kstr] != b_data[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return b_bool;
}

/* End of code generation (strcmp.c) */
