/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * str2double1.c
 *
 * Code generation for function 'str2double1'
 *
 */

/* Include files */
#include "str2double1.h"
#include "mapTF_data.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Declarations */
static void b_copysign(emxArray_char_T *s1, int *idx, const emxArray_char_T *s,
                       int *k, int n, bool *foundsign, bool *success);

static bool copydigits(emxArray_char_T *s1, int *idx, const emxArray_char_T *s,
                       int *k, int n, bool allowpoint);

static bool copyexponent(emxArray_char_T *s1, int *idx,
                         const emxArray_char_T *s, int *k, int n);

static bool isUnitImag(const emxArray_char_T *s, int k, int n);

static void readNonFinite(const emxArray_char_T *s, int *k, int n,
                          bool *b_finite, double *fv);

/* Function Definitions */
/*
 *
 */
static void b_copysign(emxArray_char_T *s1, int *idx, const emxArray_char_T *s,
                       int *k, int n, bool *foundsign, bool *success)
{
  const char *s_data;
  char c;
  char *s1_data;
  bool exitg1;
  bool isneg;
  s_data = s->data;
  s1_data = s1->data;
  isneg = false;
  *foundsign = false;
  exitg1 = false;
  while ((!exitg1) && (*k <= n)) {
    c = s_data[*k - 1];
    if (c == '-') {
      isneg = !isneg;
      *foundsign = true;
      (*k)++;
    } else if (c == ',') {
      (*k)++;
    } else if (c == '+') {
      *foundsign = true;
      (*k)++;
    } else if (!bv[(unsigned char)c & 127]) {
      exitg1 = true;
    } else {
      (*k)++;
    }
  }
  *success = (*k <= n);
  if ((*success) && isneg) {
    if ((*idx >= 2) && (s1_data[*idx - 2] == '-')) {
      s1_data[*idx - 2] = ' ';
    } else {
      s1_data[*idx - 1] = '-';
      (*idx)++;
    }
  }
}

/*
 *
 */
static bool copydigits(emxArray_char_T *s1, int *idx, const emxArray_char_T *s,
                       int *k, int n, bool allowpoint)
{
  const char *s_data;
  char c;
  char *s1_data;
  bool exitg1;
  bool haspoint;
  bool success;
  s_data = s->data;
  s1_data = s1->data;
  success = (*k <= n);
  haspoint = false;
  exitg1 = false;
  while ((!exitg1) && (success && (*k <= n))) {
    c = s_data[*k - 1];
    if ((c >= '0') && (c <= '9')) {
      s1_data[*idx - 1] = c;
      (*idx)++;
      (*k)++;
    } else if (c == '.') {
      if (allowpoint && (!haspoint)) {
        success = true;
      } else {
        success = false;
      }
      if (success) {
        s1_data[*idx - 1] = '.';
        (*idx)++;
        haspoint = true;
      }
      (*k)++;
    } else if (c == ',') {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }
  return success;
}

/*
 *
 */
static bool copyexponent(emxArray_char_T *s1, int *idx,
                         const emxArray_char_T *s, int *k, int n)
{
  int kexp;
  const char *s_data;
  char c;
  char *s1_data;
  bool b_success;
  bool success;
  s_data = s->data;
  s1_data = s1->data;
  success = true;
  if (*k <= n) {
    c = s_data[*k - 1];
    if ((c == 'E') || (c == 'e')) {
      s1_data[*idx - 1] = 'e';
      (*idx)++;
      (*k)++;
      while ((*k <= n) && (s_data[*k - 1] == ',')) {
        (*k)++;
      }
      if (*k <= n) {
        if (s_data[*k - 1] == '-') {
          s1_data[*idx - 1] = '-';
          (*idx)++;
          (*k)++;
        } else if (s_data[*k - 1] == '+') {
          (*k)++;
        }
      }
      kexp = *k;
      b_success = copydigits(s1, idx, s, k, n, false);
      if ((!b_success) || (*k <= kexp)) {
        success = false;
      }
    }
  }
  return success;
}

/*
 *
 */
static bool isUnitImag(const emxArray_char_T *s, int k, int n)
{
  int b_k;
  int j;
  char b_c[3];
  const char *s_data;
  char c;
  bool p;
  s_data = s->data;
  p = false;
  if (k <= n) {
    c = s_data[k - 1];
    if (c == 'j') {
      p = true;
    } else if (c == 'i') {
      if (k >= n - 1) {
        p = true;
      } else {
        b_k = k;
        for (j = 0; j < 3; j++) {
          b_c[j] = '\x00';
          while ((b_k <= n) && (s_data[b_k - 1] == ',')) {
            b_k++;
          }
          if (b_k <= n) {
            b_c[j] = s_data[b_k - 1];
          }
          b_k++;
        }
        if ((((!(b_c[0] == 'I')) && (!(b_c[0] == 'i'))) ||
             ((!(b_c[1] == 'N')) && (!(b_c[1] == 'n'))) ||
             ((!(b_c[2] == 'F')) && (!(b_c[2] == 'f')))) &&
            (((!(b_c[0] == 'N')) && (!(b_c[0] == 'n'))) ||
             ((!(b_c[1] == 'A')) && (!(b_c[1] == 'a'))) ||
             ((!(b_c[2] == 'N')) && (!(b_c[2] == 'n'))))) {
          p = true;
        }
      }
    }
  }
  return p;
}

/*
 *
 */
static void readNonFinite(const emxArray_char_T *s, int *k, int n,
                          bool *b_finite, double *fv)
{
  int j;
  int ksaved;
  char c[3];
  const char *s_data;
  s_data = s->data;
  ksaved = *k;
  for (j = 0; j < 3; j++) {
    c[j] = '\x00';
    while ((*k <= n) && (s_data[*k - 1] == ',')) {
      (*k)++;
    }
    if (*k <= n) {
      c[j] = s_data[*k - 1];
    }
    (*k)++;
  }
  if (((c[0] == 'I') || (c[0] == 'i')) && ((c[1] == 'N') || (c[1] == 'n')) &&
      ((c[2] == 'F') || (c[2] == 'f'))) {
    *b_finite = false;
    *fv = rtInf;
  } else if (((c[0] == 'N') || (c[0] == 'n')) &&
             ((c[1] == 'A') || (c[1] == 'a')) &&
             ((c[2] == 'N') || (c[2] == 'n'))) {
    *b_finite = false;
    *fv = rtNaN;
  } else {
    *b_finite = true;
    *fv = 0.0;
    *k = ksaved;
  }
}

/*
 *
 */
void readfloat(emxArray_char_T *s1, int *idx, const emxArray_char_T *s, int *k,
               int n, bool *isimag, bool *b_finite, double *nfv,
               bool *foundsign, bool *success)
{
  int b_idx;
  int exitg2;
  int exitg3;
  const char *s_data;
  char c;
  char *s1_data;
  bool b_foundsign;
  bool exitg4;
  s_data = s->data;
  *isimag = false;
  *b_finite = true;
  *nfv = 0.0;
  b_idx = *idx;
  b_copysign(s1, &b_idx, s, k, n, foundsign, success);
  s1_data = s1->data;
  *idx = b_idx;
  if (*success) {
    if (isUnitImag(s, *k, n)) {
      *isimag = true;
      (*k)++;
      exitg4 = false;
      while ((!exitg4) && (*k <= n)) {
        c = s_data[*k - 1];
        if (bv[(unsigned char)c & 127] || (c == '\x00') || (c == ',')) {
          (*k)++;
        } else {
          exitg4 = true;
        }
      }
      if ((*k <= n) && (s_data[*k - 1] == '*')) {
        (*k)++;
        b_copysign(s1, &b_idx, s, k, n, &b_foundsign, success);
        s1_data = s1->data;
        *idx = b_idx;
        if (*success) {
          if (isUnitImag(s, *k, n)) {
            *success = false;
          } else {
            readNonFinite(s, k, n, b_finite, nfv);
            if (*b_finite) {
              *success = copydigits(s1, idx, s, k, n, true);
              if (*success) {
                *success = copyexponent(s1, idx, s, k, n);
              }
            } else if ((b_idx >= 2) && (s1_data[b_idx - 2] == '-')) {
              *idx = b_idx - 1;
              s1_data[b_idx - 2] = ' ';
              *nfv = -*nfv;
            }
            exitg4 = false;
            while ((!exitg4) && (*k <= n)) {
              c = s_data[*k - 1];
              if (bv[(unsigned char)c & 127] || (c == '\x00') || (c == ',')) {
                (*k)++;
              } else {
                exitg4 = true;
              }
            }
            if ((*k <= n) && (s_data[*k - 1] == '*')) {
              (*k)++;
              exitg4 = false;
              while ((!exitg4) && (*k <= n)) {
                c = s_data[*k - 1];
                if (bv[(unsigned char)c & 127] || (c == '\x00') || (c == ',')) {
                  (*k)++;
                } else {
                  exitg4 = true;
                }
              }
            }
            if (*k <= n) {
              c = s_data[*k - 1];
              if ((c == 'i') || (c == 'j')) {
                (*k)++;
              }
            }
          }
          do {
            exitg3 = 0;
            if (*k <= n) {
              c = s_data[*k - 1];
              if (bv[(unsigned char)c & 127] || (c == '\x00') || (c == ',')) {
                (*k)++;
              } else {
                exitg3 = 1;
              }
            } else {
              exitg3 = 2;
            }
          } while (exitg3 == 0);
        }
      } else {
        s1_data[b_idx - 1] = '1';
        *idx = b_idx + 1;
      }
    } else {
      readNonFinite(s, k, n, b_finite, nfv);
      if (*b_finite) {
        *success = copydigits(s1, idx, s, k, n, true);
        if (*success) {
          *success = copyexponent(s1, idx, s, k, n);
        }
      } else if ((b_idx >= 2) && (s1_data[b_idx - 2] == '-')) {
        *idx = b_idx - 1;
        s1_data[b_idx - 2] = ' ';
        *nfv = -*nfv;
      }
      exitg4 = false;
      while ((!exitg4) && (*k <= n)) {
        if (bv[(unsigned char)s_data[*k - 1] & 127]) {
          (*k)++;
        } else {
          c = s_data[*k - 1];
          if ((c == '\x00') || (c == ',')) {
            (*k)++;
          } else {
            exitg4 = true;
          }
        }
      }
      if ((*k <= n) && (s_data[*k - 1] == '*')) {
        (*k)++;
        while ((*k <= n) &&
               (bv[(unsigned char)s_data[*k - 1] & 127] ||
                (s_data[*k - 1] == '\x00') || (s_data[*k - 1] == ','))) {
          (*k)++;
        }
      }
      if (*k <= n) {
        c = s_data[*k - 1];
        if ((c == 'i') || (c == 'j')) {
          (*k)++;
          *isimag = true;
        }
      }
    }
    do {
      exitg2 = 0;
      if (*k <= n) {
        if (bv[(unsigned char)s_data[*k - 1] & 127]) {
          (*k)++;
        } else {
          c = s_data[*k - 1];
          if ((c == '\x00') || (c == ',')) {
            (*k)++;
          } else {
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);
  }
}

/* End of code generation (str2double1.c) */
