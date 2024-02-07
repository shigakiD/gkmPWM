/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fileManager.c
 *
 * Code generation for function 'fileManager'
 *
 */

/* Include files */
#include "fileManager.h"
#include "getgkmweights_emxutil.h"
#include "getgkmweights_rtwutil.h"
#include "getgkmweights_types.h"
#include <stdio.h>

/* Variable Definitions */
static FILE *eml_openfiles[20];

static bool eml_autoflush[20];

/* Function Declarations */
static signed char filedata(void);

/* Function Definitions */
/*
 *
 */
static signed char filedata(void)
{
  int k;
  signed char f;
  bool exitg1;
  f = 0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 20)) {
    if (eml_openfiles[k] == NULL) {
      f = (signed char)(k + 1);
      exitg1 = true;
    } else {
      k++;
    }
  }
  return f;
}

/*
 *
 */
int cfclose(double fid)
{
  FILE *filestar;
  int cst;
  int st;
  signed char b_fileid;
  signed char fileid;
  st = -1;
  fileid = (signed char)rt_roundd(fid);
  if ((fileid < 0) || (fid != fileid)) {
    fileid = -1;
  }
  b_fileid = fileid;
  if (fileid < 0) {
    b_fileid = -1;
  }
  if (b_fileid >= 3) {
    filestar = eml_openfiles[b_fileid - 3];
  } else if (b_fileid == 0) {
    filestar = stdin;
  } else if (b_fileid == 1) {
    filestar = stdout;
  } else if (b_fileid == 2) {
    filestar = stderr;
  } else {
    filestar = NULL;
  }
  if ((filestar != NULL) && (fileid >= 3)) {
    cst = fclose(filestar);
    if (cst == 0) {
      st = 0;
      eml_openfiles[fileid - 3] = NULL;
      eml_autoflush[fileid - 3] = true;
    }
  }
  return st;
}

/*
 *
 */
signed char cfopen(const emxArray_char_T *cfilename, const char *cpermission)
{
  FILE *filestar;
  emxArray_char_T *ccfilename;
  int i;
  int loop_ub;
  const char *cfilename_data;
  signed char fileid;
  signed char j;
  char *ccfilename_data;
  cfilename_data = cfilename->data;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    emxInit_char_T(&ccfilename, 2);
    i = ccfilename->size[0] * ccfilename->size[1];
    ccfilename->size[0] = 1;
    ccfilename->size[1] = cfilename->size[1] + 1;
    emxEnsureCapacity_char_T(ccfilename, i);
    ccfilename_data = ccfilename->data;
    loop_ub = cfilename->size[1];
    for (i = 0; i < loop_ub; i++) {
      ccfilename_data[i] = cfilename_data[i];
    }
    ccfilename_data[cfilename->size[1]] = '\x00';
    filestar = fopen(&ccfilename_data[0], cpermission);
    emxFree_char_T(&ccfilename);
    if (filestar != NULL) {
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      fileid = (signed char)(j + 2);
    }
  }
  return fileid;
}

/*
 *
 */
void filedata_init(void)
{
  FILE *a;
  int i;
  a = NULL;
  for (i = 0; i < 20; i++) {
    eml_autoflush[i] = false;
    eml_openfiles[i] = a;
  }
}

/*
 *
 */
void getfilestar(double fid, FILE **filestar, bool *autoflush)
{
  signed char fileid;
  fileid = (signed char)rt_roundd(fid);
  if ((fileid < 0) || (fid != fileid)) {
    fileid = -1;
  }
  if (fileid >= 3) {
    *filestar = eml_openfiles[fileid - 3];
    *autoflush = eml_autoflush[fileid - 3];
  } else if (fileid == 0) {
    *filestar = stdin;
    *autoflush = true;
  } else if (fileid == 1) {
    *filestar = stdout;
    *autoflush = true;
  } else if (fileid == 2) {
    *filestar = stderr;
    *autoflush = true;
  } else {
    *filestar = NULL;
    *autoflush = true;
  }
}

/* End of code generation (fileManager.c) */
