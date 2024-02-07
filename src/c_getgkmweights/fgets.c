/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fgets.c
 *
 * Code generation for function 'fgets'
 *
 */

/* Include files */
#include "fgets.h"
#include "fileManager.h"
#include "getgkmweights_emxutil.h"
#include "getgkmweights_types.h"
#include "stdio.h"
#include <stddef.h>
#include <stdio.h>

/* Function Definitions */
/*
 *
 */
void b_fgets(double fileID, emxArray_char_T *line)
{
  FILE *b_NULL;
  FILE *b_filestar;
  FILE *filestar;
  char *cOut;
  int b_i;
  int carriageReturnAt;
  int exitg1;
  int i;
  int i1;
  int idx;
  int reachedEndOfFile;
  int st;
  int wherefrom;
  char ReadBuff[1024];
  unsigned char buf;
  char *line_data;
  bool exitg2;
  bool fileEndAfterCarriageReturn;
  bool newLineAfterCarriageReturn;
  bool readNewline;
  getfilestar(fileID, &filestar, &newLineAfterCarriageReturn);
  if ((fileID == 0.0) || (fileID == 1.0) || (fileID == 2.0)) {
    filestar = NULL;
  }
  line->size[0] = 1;
  line->size[1] = 0;
  if (!(filestar == NULL)) {
    do {
      exitg1 = 0;
      cOut = fgets(&ReadBuff[0], 1024, filestar);
      readNewline = false;
      b_NULL = NULL;
      getfilestar(fileID, &b_filestar, &newLineAfterCarriageReturn);
      if (b_filestar == b_NULL) {
        reachedEndOfFile = 0;
      } else {
        st = feof(b_filestar);
        reachedEndOfFile = ((int)st != 0);
      }
      if (cOut == NULL) {
        exitg1 = 1;
      } else {
        idx = 1;
        carriageReturnAt = 0;
        if (reachedEndOfFile != 0) {
          i = 0;
          exitg2 = false;
          while ((!exitg2) && (i < 1024)) {
            if (ReadBuff[i] == '\x00') {
              idx = i + 1;
              exitg2 = true;
            } else {
              if ((carriageReturnAt == 0) && (ReadBuff[i] == '\x0d')) {
                carriageReturnAt = i + 1;
              }
              i++;
            }
          }
          if (ReadBuff[idx - 1] == '\x00') {
            idx--;
          }
        } else {
          i = 0;
          exitg2 = false;
          while ((!exitg2) && (i < 1025)) {
            if (i + 1 > 1024) {
              idx = 1023;
              exitg2 = true;
            } else if (ReadBuff[i] == '\x0a') {
              idx = i + 1;
              exitg2 = true;
            } else {
              if ((carriageReturnAt == 0) && (ReadBuff[i] == '\x0d')) {
                carriageReturnAt = i + 1;
              }
              i++;
            }
          }
          readNewline = (ReadBuff[idx - 1] == '\x0a');
        }
        if ((carriageReturnAt > 0) && (carriageReturnAt < 1024)) {
          newLineAfterCarriageReturn = (ReadBuff[carriageReturnAt] == '\x0a');
          if ((reachedEndOfFile != 0) &&
              (ReadBuff[carriageReturnAt] == '\x00')) {
            fileEndAfterCarriageReturn = true;
          } else {
            fileEndAfterCarriageReturn = false;
          }
        } else {
          newLineAfterCarriageReturn = false;
          fileEndAfterCarriageReturn = false;
        }
        if ((carriageReturnAt == 0) || newLineAfterCarriageReturn ||
            fileEndAfterCarriageReturn) {
          if (1 > idx) {
            i = 0;
          } else {
            i = idx;
          }
          b_i = line->size[1];
          i1 = line->size[0] * line->size[1];
          line->size[1] += i;
          emxEnsureCapacity_char_T(line, i1);
          line_data = line->data;
          for (i1 = 0; i1 < i; i1++) {
            line_data[b_i + i1] = ReadBuff[i1];
          }
        } else {
          b_i = line->size[1];
          i1 = line->size[0] * line->size[1];
          line->size[1] += carriageReturnAt;
          emxEnsureCapacity_char_T(line, i1);
          line_data = line->data;
          for (i1 = 0; i1 < carriageReturnAt; i1++) {
            line_data[b_i + i1] = ReadBuff[i1];
          }
          wherefrom = SEEK_CUR;
          getfilestar(fileID, &b_filestar, &newLineAfterCarriageReturn);
          if ((fileID == 0.0) || (fileID == 1.0) || (fileID == 2.0)) {
            b_filestar = NULL;
          }
          if (!(b_filestar == NULL)) {
            fseek(b_filestar, (long int)(carriageReturnAt - idx), wherefrom);
          }
        }
        if (readNewline || (reachedEndOfFile != 0) || (carriageReturnAt > 0)) {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
    b_NULL = NULL;
    getfilestar(fileID, &filestar, &newLineAfterCarriageReturn);
    if (filestar == b_NULL) {
      b_i = 0;
    } else {
      st = feof(filestar);
      b_i = ((int)st != 0);
    }
    if (b_i == 0) {
      getfilestar(fileID, &filestar, &newLineAfterCarriageReturn);
      if ((fileID == 0.0) || (fileID == 1.0) || (fileID == 2.0)) {
        filestar = NULL;
      }
      if (!(filestar == NULL)) {
        fread(&buf, sizeof(unsigned char), (size_t)1, filestar);
      }
      b_NULL = NULL;
      getfilestar(fileID, &filestar, &newLineAfterCarriageReturn);
      if (filestar == b_NULL) {
        b_i = 0;
      } else {
        st = feof(filestar);
        b_i = ((int)st != 0);
      }
      if (b_i == 0) {
        wherefrom = SEEK_CUR;
        getfilestar(fileID, &filestar, &newLineAfterCarriageReturn);
        if ((fileID == 0.0) || (fileID == 1.0) || (fileID == 2.0)) {
          filestar = NULL;
        }
        if (!(filestar == NULL)) {
          fseek(filestar, (long int)-1.0, wherefrom);
        }
      }
    }
  }
}

/* End of code generation (fgets.c) */
