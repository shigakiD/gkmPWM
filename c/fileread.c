/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fileread.c
 *
 * Code generation for function 'fileread'
 *
 */

/* Include files */
#include "fileread.h"
#include "fileManager.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void fileread(const emxArray_char_T *filename, emxArray_char_T *out)
{
  static const char cv[3] = {'a', 'l', 'l'};
  FILE *filestar;
  size_t nBytes;
  size_t numReadSizeT;
  emxArray_char_T *At;
  int c;
  int exitg1;
  int fid;
  int i;
  int i1;
  int numRead;
  char tbuf[1024];
  const char *filename_data;
  signed char fileid;
  char *At_data;
  char *out_data;
  bool b_bool;
  filename_data = filename->data;
  b_bool = false;
  if (filename->size[1] == 3) {
    numRead = 0;
    do {
      exitg1 = 0;
      if (numRead < 3) {
        if (filename_data[numRead] != cv[numRead]) {
          exitg1 = 1;
        } else {
          numRead++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  if (b_bool) {
    fid = 0;
  } else {
    fileid = cfopen(filename, "rb");
    fid = fileid;
  }
  emxInit_char_T(&At, 1);
  At_data = At->data;
  nBytes = sizeof(char);
  getfilestar(fid, &filestar, &b_bool);
  if ((fid == 0) || (fid == 1) || (fid == 2)) {
    filestar = NULL;
  }
  At->size[0] = 0;
  if (!(filestar == NULL)) {
    c = 1;
    while (c > 0) {
      c = 0;
      numRead = 1;
      while ((c < 1024) && (numRead > 0)) {
        numReadSizeT = fread(&tbuf[c], nBytes, (size_t)(1024 - c), filestar);
        numRead = (int)numReadSizeT;
        c += (int)numReadSizeT;
      }
      if (1 > c) {
        numRead = 0;
      } else {
        numRead = c;
      }
      i = At->size[0];
      i1 = At->size[0];
      At->size[0] += numRead;
      emxEnsureCapacity_char_T(At, i1);
      At_data = At->data;
      for (i1 = 0; i1 < numRead; i1++) {
        At_data[i + i1] = tbuf[i1];
      }
    }
  }
  i = out->size[0] * out->size[1];
  out->size[0] = 1;
  out->size[1] = At->size[0];
  emxEnsureCapacity_char_T(out, i);
  out_data = out->data;
  numRead = At->size[0];
  for (i = 0; i < numRead; i++) {
    out_data[i] = At_data[i];
  }
  emxFree_char_T(&At);
  cfclose(fid);
}

/* End of code generation (fileread.c) */
