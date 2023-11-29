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
#include "find.h"
#include "fseek.h"
#include "ftell.h"
#include "gkmPWMlasso4_data.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include "sort.h"
#include "str2double.h"
#include "strip.h"
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 * function mat = getmotif(filename, m)
 */
void getmotif(const emxArray_char_T *filename, const emxArray_real_T *m,
              emxArray_cell_wrap_0 *mat)
{
  static const char b[5] = {'M', 'O', 'T', 'I', 'F'};
  static const char cv[3] = {'a', 'l', 'l'};
  static const char cv1[2] = {' ', '\x09'};
  cell_wrap_0 *mat_data;
  cell_wrap_0 *mat_tmp_data;
  emxArray_boolean_T *y;
  emxArray_cell_wrap_0 *mat_tmp;
  emxArray_char_T *b_fid;
  emxArray_char_T *b_x;
  emxArray_char_T *line;
  emxArray_char_T *x;
  emxArray_int32_T *b_n;
  emxArray_int32_T *ind2;
  emxArray_int32_T *r1;
  emxArray_real_T *b_r;
  emxArray_real_T *n;
  emxArray_real_T *tmp;
  creal_T dc;
  const double *m_data;
  double b_i;
  double col_count;
  double loc;
  double row_count;
  double *n_data;
  double *r3;
  double *tmp_data;
  int b_loop_ub;
  int exitg1;
  int fid;
  int i;
  int i1;
  int i2;
  int idx;
  int itoken;
  int j;
  int k;
  int loop_ub;
  int nz;
  int vlen;
  int *ind2_data;
  int *r2;
  char a[5];
  const char *filename_data;
  signed char fileid;
  char *line_data;
  char *x_data;
  bool b_bool;
  bool exitg2;
  bool *y_data;
  m_data = m->data;
  filename_data = filename->data;
  /*  filename is the meme file that contains the motifs; */
  /*  n is the nth motif in the file */
  /* 'getmotif:4' mat = cell(length(m),1); */
  vlen = m->size[1];
  i = mat->size[0];
  mat->size[0] = m->size[1];
  emxEnsureCapacity_cell_wrap_0(mat, i);
  mat_data = mat->data;
  for (i = 0; i < vlen; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 0;
  }
  emxInit_real_T(&n, 2);
  /* 'getmotif:5' mat = coder.nullcopy(mat); */
  /* 'getmotif:6' [n ind] = sort(m); */
  i = n->size[0] * n->size[1];
  n->size[0] = 1;
  n->size[1] = m->size[1];
  emxEnsureCapacity_real_T(n, i);
  n_data = n->data;
  loop_ub = m->size[1];
  for (i = 0; i < loop_ub; i++) {
    n_data[i] = m_data[i];
  }
  emxInit_int32_T(&b_n, 2);
  b_sort(n, b_n);
  n_data = n->data;
  /* 'getmotif:7' fid = fopen(filename); */
  b_bool = false;
  emxFree_int32_T(&b_n);
  if (filename->size[1] == 3) {
    vlen = 0;
    do {
      exitg1 = 0;
      if (vlen < 3) {
        if (filename_data[vlen] != cv[vlen]) {
          exitg1 = 1;
        } else {
          vlen++;
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
  /* 'getmotif:8' if fid < 0 */
  emxInit_boolean_T(&y, 2);
  emxInit_real_T(&b_r, 2);
  if (fid < 0) {
    /* 'getmotif:9' fprintf('The file cannot be opened\n') */
    printf("The file cannot be opened\n");
    fflush(stdout);
  } else {
    /* 'getmotif:10' else */
    /* 'getmotif:11' i=0; */
    b_i = 0.0;
    /* 'getmotif:12' for j = 1:length(n) */
    i = n->size[1];
    emxInit_char_T(&line, 2);
    emxInit_real_T(&tmp, 2);
    emxInit_char_T(&x, 2);
    emxInit_char_T(&b_x, 2);
    emxInit_char_T(&b_fid, 2);
    for (j = 0; j < i; j++) {
      /* 'getmotif:13' while i~=n(j) */
      while (b_i != n_data[j]) {
        /* 'getmotif:14' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
        /* 'getmotif:15' if length(line) >= 5 */
        if (line->size[1] >= 5) {
          /* 'getmotif:16' if strcmp(line(1:5), 'MOTIF') */
          for (i1 = 0; i1 < 5; i1++) {
            a[i1] = line_data[i1];
          }
          vlen = memcmp(&a[0], &b[0], 5);
          if (vlen == 0) {
            /* 'getmotif:17' i = i+1; */
            b_i++;
          }
        }
      }
      /* 'getmotif:21' line = fgetl(fid); */
      fgetl(fid, b_fid);
      /* 'getmotif:23' loc = ftell(fid); */
      loc = b_ftell(fid);
      /* 'getmotif:24' line = fgetl(fid); */
      fgetl(fid, line);
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
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = line->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        i1 = line->size[1];
        for (k = 0; k < i1; k++) {
          y_data[k] = bv[(unsigned char)line_data[k] & 127];
        }
        b_diff(y, b_r);
        r3 = b_r->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = b_r->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          y_data[i1] = (r3[i1] == 1.0);
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
        fgetl(fid, line);
        line_data = line->data;
      }
      /* 'getmotif:34' fseek(fid, loc, "bof"); */
      b_fseek(fid, loc);
      /* 'getmotif:35' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;
      /* 'getmotif:37' tmp = zeros(row_count, col_count); */
      i1 = tmp->size[0] * tmp->size[1];
      tmp->size[0] = (int)row_count;
      tmp->size[1] = (int)col_count;
      emxEnsureCapacity_real_T(tmp, i1);
      tmp_data = tmp->data;
      loop_ub = (int)row_count * (int)col_count;
      for (i1 = 0; i1 < loop_ub; i1++) {
        tmp_data[i1] = 0.0;
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
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = line->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        i1 = line->size[1];
        for (k = 0; k < i1; k++) {
          y_data[k] = bv[(unsigned char)line_data[k] & 127];
        }
        b_diff(y, b_r);
        r3 = b_r->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = b_r->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        loop_ub = b_r->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          y_data[i1] = (r3[i1] == 1.0);
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
          strip(line, x);
          x_data = x->data;
          vlen = x->size[1];
          k = 0;
          exitg2 = false;
          while ((!exitg2) && (k + 1 <= vlen)) {
            loop_ub = 0;
            do {
              exitg1 = 0;
              if (loop_ub < 2) {
                if (x_data[k] == cv1[loop_ub]) {
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
                if (x_data[k] == cv1[loop_ub]) {
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
          if (k + 1 > x->size[1]) {
            i1 = 0;
            i2 = 0;
          } else {
            i1 = k;
            i2 = x->size[1];
          }
          vlen = line->size[0] * line->size[1];
          line->size[0] = 1;
          loop_ub = i2 - i1;
          line->size[1] = loop_ub;
          emxEnsureCapacity_char_T(line, vlen);
          line_data = line->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            line_data[i2] = x_data[i1 + i2];
          }
          if (itoken > k) {
            i1 = 0;
            k = 0;
          } else {
            i1 = itoken - 1;
          }
          /* 'getmotif:47' tmp(curr_row, idx)  = real(str2double(tok)); */
          i2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          loop_ub = k - i1;
          b_x->size[1] = loop_ub;
          emxEnsureCapacity_char_T(b_x, i2);
          line_data = b_x->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            line_data[i2] = x_data[i1 + i2];
          }
          dc = str2double(b_x);
          tmp_data[((int)loc + tmp->size[0] * idx) - 1] = dc.re;
        }
        /* 'getmotif:49' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
        /* 'getmotif:50' curr_row = curr_row + 1; */
        loc++;
      }
      /* 'getmotif:52' mat{j} = tmp; */
      i1 = mat_data[j].f1->size[0] * mat_data[j].f1->size[1];
      mat_data[j].f1->size[0] = tmp->size[0];
      mat_data[j].f1->size[1] = tmp->size[1];
      emxEnsureCapacity_real_T(mat_data[j].f1, i1);
      loop_ub = tmp->size[0] * tmp->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        mat_data[j].f1->data[i1] = tmp_data[i1];
      }
    }
    emxFree_char_T(&b_fid);
    emxFree_char_T(&b_x);
    emxFree_char_T(&x);
    emxFree_real_T(&tmp);
    emxFree_char_T(&line);
  }
  /* 'getmotif:55' fclose(fid); */
  cfclose(fid);
  /* 'getmotif:56' if length(n)==1 */
  if (n->size[1] != 1) {
    emxInit_int32_T(&ind2, 1);
    /* 'getmotif:59' else */
    /* 'getmotif:60' mylength = length(n); */
    /* 'getmotif:61' ind2 = zeros(mylength,1); */
    i = ind2->size[0];
    ind2->size[0] = n->size[1];
    emxEnsureCapacity_int32_T(ind2, i);
    ind2_data = ind2->data;
    loop_ub = n->size[1];
    for (i = 0; i < loop_ub; i++) {
      ind2_data[i] = 0;
    }
    /* 'getmotif:62' for i = 1:mylength */
    i = n->size[1];
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:63' ind2(i,1) = i; */
      ind2_data[vlen] = vlen + 1;
    }
    /* 'getmotif:65' [a, b] = size(ind2); */
    /* 'getmotif:66' for i = 1:length(n) */
    i = n->size[1];
    if (0 <= n->size[1] - 1) {
      b_loop_ub = n->size[1];
    }
    emxInit_int32_T(&r1, 2);
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:67' ind2(i,1) = find(n==m(i)); */
      loc = m_data[vlen];
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = n->size[1];
      emxEnsureCapacity_boolean_T(y, i1);
      y_data = y->data;
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        y_data[i1] = (n_data[i1] == loc);
      }
      eml_find(y, r1);
      r2 = r1->data;
      i1 = b_r->size[0] * b_r->size[1];
      b_r->size[0] = 1;
      b_r->size[1] = r1->size[1];
      emxEnsureCapacity_real_T(b_r, i1);
      r3 = b_r->data;
      loop_ub = r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r3[i1] = r2[i1];
      }
      ind2_data[vlen] = (int)r3[0];
    }
    emxFree_int32_T(&r1);
    emxInit_cell_wrap_0(&mat_tmp);
    /* 'getmotif:69' mat_tmp = cell(mylength,1); */
    /* 'getmotif:70' for i = 1:mylength */
    i = n->size[1];
    i1 = mat_tmp->size[0];
    mat_tmp->size[0] = n->size[1];
    emxEnsureCapacity_cell_wrap_0(mat_tmp, i1);
    mat_tmp_data = mat_tmp->data;
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:71' mat_tmp{i} = mat{ind2(i)}; */
      i1 = mat_tmp_data[vlen].f1->size[0] * mat_tmp_data[vlen].f1->size[1];
      mat_tmp_data[vlen].f1->size[0] =
          mat_data[ind2_data[vlen] - 1].f1->size[0];
      mat_tmp_data[vlen].f1->size[1] =
          mat_data[ind2_data[vlen] - 1].f1->size[1];
      emxEnsureCapacity_real_T(mat_tmp_data[vlen].f1, i1);
      loop_ub = mat_data[ind2_data[vlen] - 1].f1->size[0] *
                mat_data[ind2_data[vlen] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        mat_tmp_data[vlen].f1->data[i1] =
            mat_data[ind2_data[vlen] - 1].f1->data[i1];
      }
    }
    emxFree_int32_T(&ind2);
    /* 'getmotif:73' [a, b] = size(mat_tmp{1}); */
    /* 'getmotif:74' for i = 1:length(n) */
    i = n->size[1];
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:75' mat{i} = mat_tmp{i}; */
      i1 = mat_data[vlen].f1->size[0] * mat_data[vlen].f1->size[1];
      mat_data[vlen].f1->size[0] = mat_tmp_data[vlen].f1->size[0];
      mat_data[vlen].f1->size[1] = mat_tmp_data[vlen].f1->size[1];
      emxEnsureCapacity_real_T(mat_data[vlen].f1, i1);
      loop_ub = mat_tmp_data[vlen].f1->size[0] * mat_tmp_data[vlen].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        mat_data[vlen].f1->data[i1] = mat_tmp_data[vlen].f1->data[i1];
      }
    }
    emxFree_cell_wrap_0(&mat_tmp);
  } else {
    /*  Alternative to cell2mat */
    /*  mat = cell2mat(mat); */
  }
  emxFree_real_T(&b_r);
  emxFree_boolean_T(&y);
  emxFree_real_T(&n);
}

/* End of code generation (getmotif.c) */
