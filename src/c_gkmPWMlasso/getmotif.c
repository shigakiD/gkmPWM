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
 * function [mat, names] = getmotif(filename, m)
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
  emxArray_char_T *b_line_clean;
  emxArray_char_T *line;
  emxArray_char_T *line_clean;
  emxArray_int32_T *b_n;
  emxArray_int32_T *ind2;
  emxArray_int32_T *match_out;
  emxArray_real_T *n;
  emxArray_real_T *tmp;
  emxArray_real_T *zero_idx;
  creal_T dc;
  const double *m_data;
  double b_i;
  double col_count;
  double loc;
  double row_count;
  double *n_data;
  double *tmp_data;
  double *zero_idx_data;
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
  int *match_out_data;
  char a[5];
  const char *filename_data;
  signed char fileid;
  char *line_clean_data;
  char *line_data;
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
  /* 'getmotif:6' names = cell(length(m),1); */
  /* 'getmotif:7' names = coder.nullcopy(names); */
  /* 'getmotif:8' [n ind] = sort(m); */
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
  d_sort(n, b_n);
  n_data = n->data;
  /* 'getmotif:9' fid = fopen(filename); */
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
  /* 'getmotif:10' if fid < 0 */
  emxInit_real_T(&zero_idx, 2);
  emxInit_boolean_T(&y, 2);
  if (fid < 0) {
    /* 'getmotif:11' fprintf('The file cannot be opened\n') */
    printf("The file cannot be opened\n");
    fflush(stdout);
  } else {
    /* 'getmotif:12' else */
    /* 'getmotif:13' i=0; */
    b_i = 0.0;
    /* 'getmotif:14' for j = 1:length(n) */
    i = n->size[1];
    emxInit_char_T(&line, 2);
    emxInit_char_T(&line_clean, 2);
    emxInit_real_T(&tmp, 2);
    emxInit_char_T(&b_line_clean, 2);
    emxInit_char_T(&b_fid, 2);
    for (j = 0; j < i; j++) {
      /* 'getmotif:15' while i~=n(j) */
      while (b_i != n_data[j]) {
        /* 'getmotif:16' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
        /* 'getmotif:17' if length(line) >= 5 */
        if (line->size[1] >= 5) {
          /* 'getmotif:18' if strcmp(line(1:5), 'MOTIF') */
          for (i1 = 0; i1 < 5; i1++) {
            a[i1] = line_data[i1];
          }
          vlen = memcmp(&a[0], &b[0], 5);
          if (vlen == 0) {
            /* 'getmotif:19' i = i+1; */
            b_i++;
            /* 'getmotif:20' line_clean = strip(line); */
            /* 'getmotif:21' zero_idx = strfind(line_clean, ' '); */
            /* 'getmotif:22' target_zero = zero_idx(end); */
            /* 'getmotif:23' names{i} = line_clean(target_zero+1:end); */
            /*  a = strsplit(line, ' '); */
            /*  names{i} = a{end}; */
          }
        }
      }
      /* 'getmotif:29' line = fgetl(fid); */
      fgetl(fid, b_fid);
      /* 'getmotif:31' loc = ftell(fid); */
      loc = b_ftell(fid);
      /* 'getmotif:32' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;
      /* 'getmotif:34' row_count = 0; */
      row_count = 0.0;
      /* 'getmotif:35' col_count = 0; */
      col_count = 0.0;
      /* 'getmotif:36' while ~isempty(line) */
      while (line->size[1] != 0) {
        /* 'getmotif:37' row_count = row_count + 1; */
        row_count++;
        /* 'getmotif:38' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:39' col_count         = sum(diff(number_whitespace)==1) +
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
        b_diff(y, zero_idx);
        zero_idx_data = zero_idx->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = zero_idx->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        loop_ub = zero_idx->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          y_data[i1] = (zero_idx_data[i1] == 1.0);
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
        /* 'getmotif:40' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
      }
      /* 'getmotif:42' fseek(fid, loc, "bof"); */
      b_fseek(fid, loc);
      /* 'getmotif:43' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;
      /* 'getmotif:45' tmp = zeros(row_count, col_count); */
      i1 = tmp->size[0] * tmp->size[1];
      tmp->size[0] = (int)row_count;
      tmp->size[1] = (int)col_count;
      emxEnsureCapacity_real_T(tmp, i1);
      tmp_data = tmp->data;
      loop_ub = (int)row_count * (int)col_count;
      for (i1 = 0; i1 < loop_ub; i1++) {
        tmp_data[i1] = 0.0;
      }
      /* 'getmotif:46' curr_row = 1; */
      loc = 1.0;
      /* 'getmotif:47' while ~isempty(line) */
      while (line->size[1] != 0) {
        /*  Alternative to str2num */
        /*  mat{j} = [mat{j}; str2num(line)]; */
        /* 'getmotif:50' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:51' number_element    = sum(diff(number_whitespace)==1) +
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
        b_diff(y, zero_idx);
        zero_idx_data = zero_idx->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = zero_idx->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        loop_ub = zero_idx->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          y_data[i1] = (zero_idx_data[i1] == 1.0);
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
        /* 'getmotif:52' remain  = line; */
        /* 'getmotif:53' for idx=1:number_element */
        for (idx = 0; idx <= nz; idx++) {
          /* 'getmotif:54' [tok, remain] = strtok(strip(remain), [char(32),
           * char(9)]); */
          strip(line, line_clean);
          line_clean_data = line_clean->data;
          vlen = line_clean->size[1];
          k = 0;
          exitg2 = false;
          while ((!exitg2) && (k + 1 <= vlen)) {
            loop_ub = 0;
            do {
              exitg1 = 0;
              if (loop_ub < 2) {
                if (line_clean_data[k] == cv1[loop_ub]) {
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
                if (line_clean_data[k] == cv1[loop_ub]) {
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
          if (k + 1 > line_clean->size[1]) {
            i1 = 0;
            i2 = 0;
          } else {
            i1 = k;
            i2 = line_clean->size[1];
          }
          vlen = line->size[0] * line->size[1];
          line->size[0] = 1;
          loop_ub = i2 - i1;
          line->size[1] = loop_ub;
          emxEnsureCapacity_char_T(line, vlen);
          line_data = line->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            line_data[i2] = line_clean_data[i1 + i2];
          }
          if (itoken > k) {
            i1 = 0;
            k = 0;
          } else {
            i1 = itoken - 1;
          }
          /* 'getmotif:55' tmp(curr_row, idx)  = real(str2double(tok)); */
          i2 = b_line_clean->size[0] * b_line_clean->size[1];
          b_line_clean->size[0] = 1;
          loop_ub = k - i1;
          b_line_clean->size[1] = loop_ub;
          emxEnsureCapacity_char_T(b_line_clean, i2);
          line_data = b_line_clean->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            line_data[i2] = line_clean_data[i1 + i2];
          }
          dc = str2double(b_line_clean);
          tmp_data[((int)loc + tmp->size[0] * idx) - 1] = dc.re;
        }
        /* 'getmotif:57' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
        /* 'getmotif:58' curr_row = curr_row + 1; */
        loc++;
      }
      /* 'getmotif:60' mat{j} = tmp; */
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
    emxFree_char_T(&b_line_clean);
    emxFree_real_T(&tmp);
    emxFree_char_T(&line_clean);
    emxFree_char_T(&line);
  }
  /* 'getmotif:63' fclose(fid); */
  cfclose(fid);
  /* 'getmotif:64' if length(n)==1 */
  if (n->size[1] != 1) {
    emxInit_int32_T(&ind2, 1);
    /* 'getmotif:67' else */
    /* 'getmotif:68' mylength = length(n); */
    /* 'getmotif:69' ind2 = zeros(mylength,1); */
    i = ind2->size[0];
    ind2->size[0] = n->size[1];
    emxEnsureCapacity_int32_T(ind2, i);
    ind2_data = ind2->data;
    loop_ub = n->size[1];
    for (i = 0; i < loop_ub; i++) {
      ind2_data[i] = 0;
    }
    /* 'getmotif:70' for i = 1:mylength */
    i = n->size[1];
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:71' ind2(i,1) = i; */
      ind2_data[vlen] = vlen + 1;
    }
    /* 'getmotif:73' [a, b] = size(ind2); */
    /* 'getmotif:74' for i = 1:length(n) */
    i = n->size[1];
    if (0 <= n->size[1] - 1) {
      b_loop_ub = n->size[1];
    }
    emxInit_int32_T(&match_out, 2);
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:75' ind2(i,1) = find(n==m(i)); */
      loc = m_data[vlen];
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = n->size[1];
      emxEnsureCapacity_boolean_T(y, i1);
      y_data = y->data;
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        y_data[i1] = (n_data[i1] == loc);
      }
      eml_find(y, match_out);
      match_out_data = match_out->data;
      i1 = zero_idx->size[0] * zero_idx->size[1];
      zero_idx->size[0] = 1;
      zero_idx->size[1] = match_out->size[1];
      emxEnsureCapacity_real_T(zero_idx, i1);
      zero_idx_data = zero_idx->data;
      loop_ub = match_out->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        zero_idx_data[i1] = match_out_data[i1];
      }
      ind2_data[vlen] = (int)zero_idx_data[0];
    }
    emxFree_int32_T(&match_out);
    emxInit_cell_wrap_0(&mat_tmp);
    /* 'getmotif:77' mat_tmp = cell(mylength,1); */
    /* 'getmotif:78' for i = 1:mylength */
    i = n->size[1];
    i1 = mat_tmp->size[0];
    mat_tmp->size[0] = n->size[1];
    emxEnsureCapacity_cell_wrap_0(mat_tmp, i1);
    mat_tmp_data = mat_tmp->data;
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:79' mat_tmp{i} = mat{ind2(i)}; */
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
    /* 'getmotif:81' [a, b] = size(mat_tmp{1}); */
    /* 'getmotif:82' for i = 1:length(n) */
    i = n->size[1];
    for (vlen = 0; vlen < i; vlen++) {
      /* 'getmotif:83' mat{i} = mat_tmp{i}; */
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
  emxFree_boolean_T(&y);
  emxFree_real_T(&n);
  emxFree_real_T(&zero_idx);
}

/* End of code generation (getmotif.c) */
