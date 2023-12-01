/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getmotif.c
 *
 * Code generation for function 'getmotif'
 *
 */

/* Include files */
#include "getmotif.h"
#include "diff.h"
#include "fgetl.h"
#include "fileManager.h"
#include "gkmPWM_data.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "sort.h"
#include "str2double.h"
#include "strip.h"
#include "lapacke.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 * function mat = getmotif(filename, m)
 */
void getmotif(cell_wrap_3 mat[1968])
{
  static const char b[5] = {'M', 'O', 'T', 'I', 'F'};
  static const char cv[2] = {' ', '\x09'};
  FILE *filestar;
  int wherefrom;
  long position_t;
  cell_wrap_3 mat_tmp[1968];
  emxArray_boolean_T *y;
  emxArray_char_T *b_fileid;
  emxArray_char_T *b_x;
  emxArray_char_T *c_x;
  emxArray_char_T *line;
  emxArray_real_T *r;
  emxArray_real_T *tmp;
  creal_T dc;
  double x[1968];
  double b_i;
  double col_count;
  double loc;
  double row_count;
  double *r1;
  double *tmp_data;
  int unusedExpr[1968];
  int exitg1;
  int i;
  int i1;
  int idx;
  int itoken;
  int j;
  int k;
  int loop_ub;
  int nz;
  int vlen;
  short ii_data[1968];
  short ind2[1968];
  char b_a[5];
  signed char fileid;
  char *line_data;
  char *x_data;
  bool d_x[1968];
  bool a;
  bool exitg2;
  bool *y_data;
  /*  filename is the meme file that contains the motifs; */
  /*  n is the nth motif in the file */
  /* 'getmotif:4' mat = cell(length(m),1); */
  /* 'getmotif:5' mat = coder.nullcopy(mat); */
  /* 'getmotif:6' [n ind] = sort(m); */
  for (i = 0; i < 1968; i++) {
    x[i] = (double)i + 1.0;
  }
  e_sort(x, unusedExpr);
  /* 'getmotif:7' fid = fopen(filename); */
  fileid = b_cfopen("combined_db_v4.meme");
  /* 'getmotif:8' if fid < 0 */
  if (fileid < 0) {
    /* 'getmotif:9' fprintf('The file cannot be opened\n') */
    printf("The file cannot be opened\n");
    fflush(stdout);
  } else {
    /* 'getmotif:10' else */
    /* 'getmotif:11' i=0; */
    b_i = 0.0;
    /* 'getmotif:12' for j = 1:length(n) */
    emxInit_char_T(&line, 2);
    emxInit_real_T(&tmp, 2);
    emxInit_char_T(&b_x, 2);
    emxInit_boolean_T(&y, 2);
    emxInit_real_T(&r, 2);
    emxInit_char_T(&c_x, 2);
    emxInit_char_T(&b_fileid, 2);
    wherefrom = SEEK_SET;
    for (j = 0; j < 1968; j++) {
      /* 'getmotif:13' while i~=n(j) */
      while (b_i != x[j]) {
        /* 'getmotif:14' line = fgetl(fid); */
        fgetl(fileid, line);
        line_data = line->data;
        /* 'getmotif:15' if length(line) >= 5 */
        if (line->size[1] >= 5) {
          /* 'getmotif:16' if strcmp(line(1:5), 'MOTIF') */
          for (i = 0; i < 5; i++) {
            b_a[i] = line_data[i];
          }
          vlen = memcmp(&b_a[0], &b[0], 5);
          if (vlen == 0) {
            /* 'getmotif:17' i = i+1; */
            b_i++;
          }
        }
      }
      /* 'getmotif:21' line = fgetl(fid); */
      fgetl(fileid, b_fileid);
      /* 'getmotif:23' loc = ftell(fid); */
      getfilestar(fileid, &filestar, &a);
      if ((fileid == 0) || (fileid == 1) || (fileid == 2)) {
        filestar = NULL;
      }
      if (filestar == NULL) {
        loc = -1.0;
      } else {
        position_t = ftell(filestar);
        loc = (double)position_t;
      }
      /* 'getmotif:24' line = fgetl(fid); */
      fgetl(fileid, line);
      line_data = line->data;
      /* 'getmotif:26' row_count = 0; */
      row_count = 0.0;
      /* 'getmotif:27' col_count = 0; */
      col_count = 0.0;
      /* 'getmotif:28' while ~isempty(line) */
      while (line->size[1] != 0) {
        /* 'getmotif:29' row_count = row_count + 1; */
        row_count++;
        /* 'getmotif:30' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:31' col_count         = sum(diff(number_whitespace)==1) +
         * 1; */
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = line->size[1];
        emxEnsureCapacity_boolean_T(y, i);
        y_data = y->data;
        i = line->size[1];
        for (k = 0; k < i; k++) {
          y_data[k] = bv[(unsigned char)line_data[k] & 127];
        }
        d_diff(y, r);
        r1 = r->data;
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = r->size[1];
        emxEnsureCapacity_boolean_T(y, i);
        y_data = y->data;
        loop_ub = r->size[1];
        for (i = 0; i < loop_ub; i++) {
          y_data[i] = (r1[i] == 1.0);
        }
        vlen = y->size[1];
        if (y->size[1] == 0) {
          nz = 0;
        } else {
          nz = y_data[0];
          for (k = 2; k <= vlen; k++) {
            nz += y_data[k - 1];
          }
        }
        col_count = (double)nz + 1.0;
        /* 'getmotif:32' line = fgetl(fid); */
        fgetl(fileid, line);
        line_data = line->data;
      }
      /* 'getmotif:34' fseek(fid, loc, "bof"); */
      if (floor(loc) == loc) {
        getfilestar(fileid, &filestar, &a);
        if ((fileid == 0) || (fileid == 1) || (fileid == 2)) {
          filestar = NULL;
        }
        if (!(filestar == NULL)) {
          fseek(filestar, (long int)loc, wherefrom);
        }
      }
      /* 'getmotif:35' line = fgetl(fid); */
      fgetl(fileid, line);
      line_data = line->data;
      /* 'getmotif:37' tmp = zeros(row_count, col_count); */
      i = tmp->size[0] * tmp->size[1];
      tmp->size[0] = (int)row_count;
      tmp->size[1] = (int)col_count;
      emxEnsureCapacity_real_T(tmp, i);
      tmp_data = tmp->data;
      loop_ub = (int)row_count * (int)col_count;
      for (i = 0; i < loop_ub; i++) {
        tmp_data[i] = 0.0;
      }
      /* 'getmotif:38' curr_row = 1; */
      loc = 1.0;
      /* 'getmotif:39' while ~isempty(line) */
      while (line->size[1] != 0) {
        /*  Alternative to str2num */
        /*  mat{j} = [mat{j}; str2num(line)]; */
        /* 'getmotif:42' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:43' number_element    = sum(diff(number_whitespace)==1) +
         * 1; */
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = line->size[1];
        emxEnsureCapacity_boolean_T(y, i);
        y_data = y->data;
        i = line->size[1];
        for (k = 0; k < i; k++) {
          y_data[k] = bv[(unsigned char)line_data[k] & 127];
        }
        d_diff(y, r);
        r1 = r->data;
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = r->size[1];
        emxEnsureCapacity_boolean_T(y, i);
        y_data = y->data;
        loop_ub = r->size[1];
        for (i = 0; i < loop_ub; i++) {
          y_data[i] = (r1[i] == 1.0);
        }
        vlen = y->size[1];
        if (y->size[1] == 0) {
          nz = 0;
        } else {
          nz = y_data[0];
          for (k = 2; k <= vlen; k++) {
            nz += y_data[k - 1];
          }
        }
        /* 'getmotif:44' remain  = line; */
        /* 'getmotif:45' for idx=1:number_element */
        for (idx = 0; idx <= nz; idx++) {
          /* 'getmotif:46' [tok, remain] = strtok(strip(remain), [char(32),
           * char(9)]); */
          strip(line, b_x);
          x_data = b_x->data;
          vlen = b_x->size[1];
          k = 0;
          exitg2 = false;
          while ((!exitg2) && (k + 1 <= vlen)) {
            loop_ub = 0;
            do {
              exitg1 = 0;
              if (loop_ub < 2) {
                if (x_data[k] == cv[loop_ub]) {
                  k++;
                  exitg1 = 1;
                } else {
                  loop_ub++;
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
          while ((!exitg2) && (k + 1 <= vlen)) {
            loop_ub = 0;
            do {
              exitg1 = 0;
              if (loop_ub < 2) {
                if (x_data[k] == cv[loop_ub]) {
                  exitg1 = 1;
                } else {
                  loop_ub++;
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
          if (k + 1 > b_x->size[1]) {
            i = 0;
            i1 = 0;
          } else {
            i = k;
            i1 = b_x->size[1];
          }
          vlen = line->size[0] * line->size[1];
          line->size[0] = 1;
          loop_ub = i1 - i;
          line->size[1] = loop_ub;
          emxEnsureCapacity_char_T(line, vlen);
          line_data = line->data;
          for (i1 = 0; i1 < loop_ub; i1++) {
            line_data[i1] = x_data[i + i1];
          }
          if (itoken > k) {
            i = 0;
            k = 0;
          } else {
            i = itoken - 1;
          }
          /* 'getmotif:47' tmp(curr_row, idx)  = real(str2double(tok)); */
          i1 = c_x->size[0] * c_x->size[1];
          c_x->size[0] = 1;
          loop_ub = k - i;
          c_x->size[1] = loop_ub;
          emxEnsureCapacity_char_T(c_x, i1);
          line_data = c_x->data;
          for (i1 = 0; i1 < loop_ub; i1++) {
            line_data[i1] = x_data[i + i1];
          }
          dc = str2double(c_x);
          tmp_data[((int)loc + tmp->size[0] * idx) - 1] = dc.re;
        }
        /* 'getmotif:49' line = fgetl(fid); */
        fgetl(fileid, line);
        line_data = line->data;
        /* 'getmotif:50' curr_row = curr_row + 1; */
        loc++;
      }
      /* 'getmotif:52' mat{j} = tmp; */
      i = mat[j].f1->size[0] * mat[j].f1->size[1];
      mat[j].f1->size[0] = tmp->size[0];
      mat[j].f1->size[1] = tmp->size[1];
      emxEnsureCapacity_real_T(mat[j].f1, i);
      loop_ub = tmp->size[0] * tmp->size[1];
      for (i = 0; i < loop_ub; i++) {
        mat[j].f1->data[i] = tmp_data[i];
      }
    }
    emxFree_char_T(&b_fileid);
    emxFree_char_T(&c_x);
    emxFree_real_T(&r);
    emxFree_boolean_T(&y);
    emxFree_char_T(&b_x);
    emxFree_real_T(&tmp);
    emxFree_char_T(&line);
  }
  emxInitMatrix_cell_wrap_3(mat_tmp);
  /* 'getmotif:55' fclose(fid); */
  cfclose(fileid);
  /* 'getmotif:56' if length(n)==1 */
  /* 'getmotif:59' else */
  /* 'getmotif:60' mylength = length(n); */
  /* 'getmotif:61' ind2 = zeros(mylength,1); */
  /* 'getmotif:62' for i = 1:mylength */
  /* 'getmotif:65' [a, b] = size(ind2); */
  /* 'getmotif:66' for i = 1:length(n) */
  /* 'getmotif:69' mat_tmp = cell(mylength,1); */
  /* 'getmotif:70' for i = 1:mylength */
  for (itoken = 0; itoken < 1968; itoken++) {
    /* 'getmotif:67' ind2(i,1) = find(n==m(i)); */
    for (i = 0; i < 1968; i++) {
      d_x[i] = (x[i] == (double)itoken + 1.0);
    }
    idx = 0;
    vlen = 0;
    exitg2 = false;
    while ((!exitg2) && (vlen < 1968)) {
      if (d_x[vlen]) {
        idx++;
        ii_data[idx - 1] = (short)(vlen + 1);
        if (idx >= 1968) {
          exitg2 = true;
        } else {
          vlen++;
        }
      } else {
        vlen++;
      }
    }
    ind2[itoken] = ii_data[0];
    /* 'getmotif:71' mat_tmp{i} = mat{ind2(i)}; */
    i = mat_tmp[itoken].f1->size[0] * mat_tmp[itoken].f1->size[1];
    i1 = ind2[itoken] - 1;
    mat_tmp[itoken].f1->size[0] = mat[i1].f1->size[0];
    mat_tmp[itoken].f1->size[1] = mat[i1].f1->size[1];
    emxEnsureCapacity_real_T(mat_tmp[itoken].f1, i);
    loop_ub = mat[i1].f1->size[0] * mat[i1].f1->size[1];
    for (i = 0; i < loop_ub; i++) {
      mat_tmp[itoken].f1->data[i] = mat[i1].f1->data[i];
    }
  }
  /* 'getmotif:73' [a, b] = size(mat_tmp{1}); */
  /* 'getmotif:74' for i = 1:length(n) */
  for (itoken = 0; itoken < 1968; itoken++) {
    /* 'getmotif:75' mat{i} = mat_tmp{i}; */
    i = mat[itoken].f1->size[0] * mat[itoken].f1->size[1];
    mat[itoken].f1->size[0] = mat_tmp[itoken].f1->size[0];
    mat[itoken].f1->size[1] = mat_tmp[itoken].f1->size[1];
    emxEnsureCapacity_real_T(mat[itoken].f1, i);
    loop_ub = mat_tmp[itoken].f1->size[0] * mat_tmp[itoken].f1->size[1];
    for (i = 0; i < loop_ub; i++) {
      mat[itoken].f1->data[i] = mat_tmp[itoken].f1->data[i];
    }
  }
  emxFreeMatrix_cell_wrap_3(mat_tmp);
}

/* End of code generation (getmotif.c) */
