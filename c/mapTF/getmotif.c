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
#include "mapTF_data.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "str2double.h"
#include "strip.h"
#include "strtok.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Function Definitions */
/*
 * function [mat, names] = getmotif(filename, m)
 */
void b_getmotif(const emxArray_char_T *filename, const emxArray_real_T *m,
                emxArray_cell_wrap_6 *mat, emxArray_cell_wrap_4 *names)
{
  static const char b[5] = {'M', 'O', 'T', 'I', 'F'};
  static const char b_cv[3] = {'a', 'l', 'l'};
  cell_wrap_4 *names_data;
  cell_wrap_6 *mat_data;
  cell_wrap_6 *mat_tmp_data;
  emxArray_boolean_T *y;
  emxArray_cell_wrap_6 *mat_tmp;
  emxArray_char_T *b_fid;
  emxArray_char_T *line;
  emxArray_char_T *line_clean;
  emxArray_char_T *tok;
  emxArray_int32_T *b_n;
  emxArray_int32_T *ind2;
  emxArray_int32_T *match_out;
  emxArray_int32_T *matches;
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
  int c_i;
  int exitg1;
  int fid;
  int i;
  int i1;
  int j;
  int loop_ub;
  int match_idx;
  int text_len;
  int *match_out_data;
  int *matches_data;
  char a[5];
  const char *filename_data;
  signed char fileid;
  char *line_clean_data;
  bool b_bool;
  bool *y_data;
  m_data = m->data;
  filename_data = filename->data;
  /*  filename is the meme file that contains the motifs; */
  /*  n is the nth motif in the file */
  /* 'getmotif:4' mat = cell(length(m),1); */
  text_len = m->size[1];
  i = mat->size[0];
  mat->size[0] = m->size[1];
  emxEnsureCapacity_cell_wrap_6(mat, i);
  mat_data = mat->data;
  for (i = 0; i < text_len; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 0;
  }
  /* 'getmotif:5' mat = coder.nullcopy(mat); */
  /* 'getmotif:6' names = cell(length(m),1); */
  text_len = m->size[1];
  i = names->size[0];
  names->size[0] = m->size[1];
  emxEnsureCapacity_cell_wrap_4(names, i);
  names_data = names->data;
  for (i = 0; i < text_len; i++) {
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }
  emxInit_real_T(&n, 2);
  /* 'getmotif:7' names = coder.nullcopy(names); */
  /* 'getmotif:8' [n ind] = sort(m); */
  i = n->size[0] * n->size[1];
  n->size[0] = 1;
  n->size[1] = m->size[1];
  emxEnsureCapacity_real_T(n, i);
  n_data = n->data;
  match_idx = m->size[1];
  for (i = 0; i < match_idx; i++) {
    n_data[i] = m_data[i];
  }
  emxInit_int32_T(&b_n, 2);
  b_sort(n, b_n);
  n_data = n->data;
  /* 'getmotif:9' fid = fopen(filename); */
  b_bool = false;
  emxFree_int32_T(&b_n);
  if (filename->size[1] == 3) {
    text_len = 0;
    do {
      exitg1 = 0;
      if (text_len < 3) {
        if (filename_data[text_len] != b_cv[text_len]) {
          exitg1 = 1;
        } else {
          text_len++;
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
  emxInit_char_T(&line, 2);
  /* 'getmotif:10' line = ''; */
  line->size[0] = 1;
  line->size[1] = 0;
  /* 'getmotif:11' if fid < 0 */
  emxInit_real_T(&zero_idx, 2);
  zero_idx_data = zero_idx->data;
  emxInit_int32_T(&match_out, 2);
  emxInit_boolean_T(&y, 2);
  if (fid < 0) {
    /* 'getmotif:12' fprintf("ERROR: Cannot open motif database.\n"); */
    printf("ERROR: Cannot open motif database.\n");
    fflush(stdout);
  } else {
    /* 'getmotif:13' else */
    /* 'getmotif:14' i=0; */
    b_i = 0.0;
    /* 'getmotif:15' for j = 1:length(n) */
    i = n->size[1];
    emxInit_char_T(&line_clean, 2);
    emxInit_real_T(&tmp, 2);
    emxInit_char_T(&tok, 2);
    emxInit_int32_T(&matches, 2);
    emxInit_char_T(&b_fid, 2);
    for (j = 0; j < i; j++) {
      /* 'getmotif:16' while i~=n(j) */
      while (b_i != n_data[j]) {
        /* 'getmotif:17' line = fgetl(fid); */
        fgetl(fid, line);
        line_clean_data = line->data;
        /* 'getmotif:18' if length(line) >= 5 */
        if (line->size[1] >= 5) {
          /* 'getmotif:19' if strcmp(line(1:5), 'MOTIF') */
          for (i1 = 0; i1 < 5; i1++) {
            a[i1] = line_clean_data[i1];
          }
          text_len = memcmp(&a[0], &b[0], 5);
          if (text_len == 0) {
            /* 'getmotif:20' i = i+1; */
            b_i++;
          }
        }
      }
      /* 'getmotif:24' line_clean = strip(line); */
      strip(line, line_clean);
      line_clean_data = line_clean->data;
      /* 'getmotif:25' zero_idx = strfind(line_clean, ' '); */
      if (line_clean->size[1] == 0) {
        zero_idx->size[0] = 1;
        zero_idx->size[1] = 0;
      } else {
        text_len = line_clean->size[1];
        i1 = matches->size[0] * matches->size[1];
        matches->size[0] = 1;
        matches->size[1] = line_clean->size[1];
        emxEnsureCapacity_int32_T(matches, i1);
        matches_data = matches->data;
        match_idx = 0;
        for (c_i = 0; c_i < text_len; c_i++) {
          if (line_clean_data[c_i] == ' ') {
            matches_data[match_idx] = c_i + 1;
            match_idx++;
          }
        }
        i1 = match_out->size[0] * match_out->size[1];
        match_out->size[0] = 1;
        match_out->size[1] = match_idx;
        emxEnsureCapacity_int32_T(match_out, i1);
        match_out_data = match_out->data;
        for (c_i = 0; c_i < match_idx; c_i++) {
          match_out_data[c_i] = matches_data[c_i];
        }
        i1 = zero_idx->size[0] * zero_idx->size[1];
        zero_idx->size[0] = 1;
        zero_idx->size[1] = match_out->size[1];
        emxEnsureCapacity_real_T(zero_idx, i1);
        zero_idx_data = zero_idx->data;
        match_idx = match_out->size[1];
        for (i1 = 0; i1 < match_idx; i1++) {
          zero_idx_data[i1] = match_out_data[i1];
        }
      }
      /* 'getmotif:26' target_zero = zero_idx(end); */
      /* 'getmotif:27' names{j} = line_clean(target_zero+1:end); */
      if ((unsigned int)zero_idx_data[zero_idx->size[1] - 1] + 1U >
          (unsigned int)line_clean->size[1]) {
        i1 = 0;
        c_i = 0;
      } else {
        i1 = (int)(unsigned int)zero_idx_data[zero_idx->size[1] - 1];
        c_i = line_clean->size[1];
      }
      text_len = names_data[j].f1->size[0] * names_data[j].f1->size[1];
      names_data[j].f1->size[0] = 1;
      match_idx = c_i - i1;
      names_data[j].f1->size[1] = match_idx;
      emxEnsureCapacity_char_T(names_data[j].f1, text_len);
      for (c_i = 0; c_i < match_idx; c_i++) {
        names_data[j].f1->data[c_i] = line_clean_data[i1 + c_i];
      }
      /*  a = strsplit(line, ' '); */
      /*  names{i} = a{end}; */
      /* 'getmotif:30' line = fgetl(fid); */
      fgetl(fid, b_fid);
      /* 'getmotif:32' loc = ftell(fid); */
      loc = b_ftell(fid);
      /* 'getmotif:33' line = fgetl(fid); */
      fgetl(fid, line);
      line_clean_data = line->data;
      /* 'getmotif:35' row_count = 0; */
      row_count = 0.0;
      /* 'getmotif:36' col_count = 0; */
      col_count = 0.0;
      /* 'getmotif:37' while ~isempty(line) */
      while (line->size[1] != 0) {
        /* 'getmotif:38' row_count = row_count + 1; */
        row_count++;
        /* 'getmotif:39' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:40' col_count         = sum(diff(number_whitespace)==1) +
         * 1; */
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = line->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        i1 = line->size[1];
        for (c_i = 0; c_i < i1; c_i++) {
          y_data[c_i] = bv[(unsigned char)line_clean_data[c_i] & 127];
        }
        diff(y, zero_idx);
        zero_idx_data = zero_idx->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = zero_idx->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        match_idx = zero_idx->size[1];
        for (i1 = 0; i1 < match_idx; i1++) {
          y_data[i1] = (zero_idx_data[i1] == 1.0);
        }
        text_len = y->size[1];
        if (y->size[1] == 0) {
          match_idx = 0;
        } else {
          match_idx = y_data[0];
          for (c_i = 2; c_i <= text_len; c_i++) {
            match_idx += y_data[c_i - 1];
          }
        }
        col_count = (double)match_idx + 1.0;
        /* 'getmotif:41' line = fgetl(fid); */
        fgetl(fid, line);
        line_clean_data = line->data;
      }
      /* 'getmotif:43' fseek(fid, loc, "bof"); */
      b_fseek(fid, loc);
      /* 'getmotif:44' line = fgetl(fid); */
      fgetl(fid, line);
      line_clean_data = line->data;
      /* 'getmotif:46' tmp = zeros(row_count, col_count); */
      i1 = tmp->size[0] * tmp->size[1];
      tmp->size[0] = (int)row_count;
      tmp->size[1] = (int)col_count;
      emxEnsureCapacity_real_T(tmp, i1);
      tmp_data = tmp->data;
      match_idx = (int)row_count * (int)col_count;
      for (i1 = 0; i1 < match_idx; i1++) {
        tmp_data[i1] = 0.0;
      }
      /* 'getmotif:47' curr_row = 1; */
      loc = 1.0;
      /* 'getmotif:48' while ~isempty(line) */
      while (line->size[1] != 0) {
        /*  Alternative to str2num */
        /*  mat{j} = [mat{j}; str2num(line)]; */
        /* 'getmotif:51' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:52' number_element    = sum(diff(number_whitespace)==1) +
         * 1; */
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = line->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        i1 = line->size[1];
        for (c_i = 0; c_i < i1; c_i++) {
          y_data[c_i] = bv[(unsigned char)line_clean_data[c_i] & 127];
        }
        diff(y, zero_idx);
        zero_idx_data = zero_idx->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = zero_idx->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        match_idx = zero_idx->size[1];
        for (i1 = 0; i1 < match_idx; i1++) {
          y_data[i1] = (zero_idx_data[i1] == 1.0);
        }
        text_len = y->size[1];
        if (y->size[1] == 0) {
          match_idx = 0;
        } else {
          match_idx = y_data[0];
          for (c_i = 2; c_i <= text_len; c_i++) {
            match_idx += y_data[c_i - 1];
          }
        }
        /* 'getmotif:53' remain  = line; */
        /* 'getmotif:54' for idx=1:number_element */
        for (text_len = 0; text_len <= match_idx; text_len++) {
          /* 'getmotif:55' [tok, remain] = strtok(strip(remain), [char(32),
           * char(9)]); */
          strip(line, line_clean);
          d_strtok(line_clean, tok, line);
          /* 'getmotif:56' tmp(curr_row, idx)  = real(str2double(tok)); */
          dc = str2double(tok);
          tmp_data[((int)loc + tmp->size[0] * text_len) - 1] = dc.re;
        }
        /* 'getmotif:58' line = fgetl(fid); */
        fgetl(fid, line);
        line_clean_data = line->data;
        /* 'getmotif:59' curr_row = curr_row + 1; */
        loc++;
      }
      /* 'getmotif:61' mat{j} = tmp; */
      i1 = mat_data[j].f1->size[0] * mat_data[j].f1->size[1];
      mat_data[j].f1->size[0] = tmp->size[0];
      mat_data[j].f1->size[1] = tmp->size[1];
      emxEnsureCapacity_real_T(mat_data[j].f1, i1);
      match_idx = tmp->size[0] * tmp->size[1];
      for (i1 = 0; i1 < match_idx; i1++) {
        mat_data[j].f1->data[i1] = tmp_data[i1];
      }
    }
    emxFree_char_T(&b_fid);
    emxFree_int32_T(&matches);
    emxFree_char_T(&tok);
    emxFree_real_T(&tmp);
    emxFree_char_T(&line_clean);
  }
  emxFree_char_T(&line);
  /* 'getmotif:64' fclose(fid); */
  cfclose(fid);
  /* 'getmotif:65' if length(n)==1 */
  if (n->size[1] != 1) {
    emxInit_int32_T(&ind2, 1);
    /* 'getmotif:68' else */
    /* 'getmotif:69' mylength = length(n); */
    /* 'getmotif:70' ind2 = zeros(mylength,1); */
    i = ind2->size[0];
    ind2->size[0] = n->size[1];
    emxEnsureCapacity_int32_T(ind2, i);
    matches_data = ind2->data;
    match_idx = n->size[1];
    for (i = 0; i < match_idx; i++) {
      matches_data[i] = 0;
    }
    /* 'getmotif:71' for i = 1:mylength */
    i = n->size[1];
    for (c_i = 0; c_i < i; c_i++) {
      /* 'getmotif:72' ind2(i,1) = i; */
      matches_data[c_i] = c_i + 1;
    }
    /* 'getmotif:74' [a, b] = size(ind2); */
    /* 'getmotif:75' for i = 1:length(n) */
    i = n->size[1];
    if (0 <= n->size[1] - 1) {
      loop_ub = n->size[1];
    }
    for (c_i = 0; c_i < i; c_i++) {
      /* 'getmotif:76' ind2(i,1) = find(n==m(i)); */
      loc = m_data[c_i];
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = n->size[1];
      emxEnsureCapacity_boolean_T(y, i1);
      y_data = y->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        y_data[i1] = (n_data[i1] == loc);
      }
      b_eml_find(y, match_out);
      match_out_data = match_out->data;
      i1 = zero_idx->size[0] * zero_idx->size[1];
      zero_idx->size[0] = 1;
      zero_idx->size[1] = match_out->size[1];
      emxEnsureCapacity_real_T(zero_idx, i1);
      zero_idx_data = zero_idx->data;
      match_idx = match_out->size[1];
      for (i1 = 0; i1 < match_idx; i1++) {
        zero_idx_data[i1] = match_out_data[i1];
      }
      matches_data[c_i] = (int)zero_idx_data[0];
    }
    emxInit_cell_wrap_6(&mat_tmp, 1);
    /* 'getmotif:78' mat_tmp = cell(mylength,1); */
    /* 'getmotif:79' for i = 1:mylength */
    i = n->size[1];
    i1 = mat_tmp->size[0];
    mat_tmp->size[0] = n->size[1];
    emxEnsureCapacity_cell_wrap_6(mat_tmp, i1);
    mat_tmp_data = mat_tmp->data;
    for (c_i = 0; c_i < i; c_i++) {
      /* 'getmotif:80' mat_tmp{i} = mat{ind2(i)}; */
      i1 = mat_tmp_data[c_i].f1->size[0] * mat_tmp_data[c_i].f1->size[1];
      mat_tmp_data[c_i].f1->size[0] =
          mat_data[matches_data[c_i] - 1].f1->size[0];
      mat_tmp_data[c_i].f1->size[1] =
          mat_data[matches_data[c_i] - 1].f1->size[1];
      emxEnsureCapacity_real_T(mat_tmp_data[c_i].f1, i1);
      match_idx = mat_data[matches_data[c_i] - 1].f1->size[0] *
                  mat_data[matches_data[c_i] - 1].f1->size[1];
      for (i1 = 0; i1 < match_idx; i1++) {
        mat_tmp_data[c_i].f1->data[i1] =
            mat_data[matches_data[c_i] - 1].f1->data[i1];
      }
    }
    emxFree_int32_T(&ind2);
    /* 'getmotif:82' [a, b] = size(mat_tmp{1}); */
    /* 'getmotif:83' for i = 1:length(n) */
    i = n->size[1];
    for (c_i = 0; c_i < i; c_i++) {
      /* 'getmotif:84' mat{i} = mat_tmp{i}; */
      i1 = mat_data[c_i].f1->size[0] * mat_data[c_i].f1->size[1];
      mat_data[c_i].f1->size[0] = mat_tmp_data[c_i].f1->size[0];
      mat_data[c_i].f1->size[1] = mat_tmp_data[c_i].f1->size[1];
      emxEnsureCapacity_real_T(mat_data[c_i].f1, i1);
      match_idx = mat_tmp_data[c_i].f1->size[0] * mat_tmp_data[c_i].f1->size[1];
      for (i1 = 0; i1 < match_idx; i1++) {
        mat_data[c_i].f1->data[i1] = mat_tmp_data[c_i].f1->data[i1];
      }
    }
    emxFree_cell_wrap_6(&mat_tmp);
  } else {
    /*  Alternative to cell2mat */
    /*  mat = cell2mat(mat); */
  }
  emxFree_boolean_T(&y);
  emxFree_int32_T(&match_out);
  emxFree_real_T(&n);
  emxFree_real_T(&zero_idx);
}

/*
 * function [mat, names] = getmotif(filename, m)
 */
void getmotif(const emxArray_char_T *filename, const emxArray_real_T *m,
              emxArray_cell_wrap_6 *mat)
{
  static const char b[5] = {'M', 'O', 'T', 'I', 'F'};
  static const char b_cv[3] = {'a', 'l', 'l'};
  FILE *filestar;
  long position_t;
  cell_wrap_6 *mat_data;
  cell_wrap_6 *mat_tmp_data;
  emxArray_boolean_T *c_n;
  emxArray_boolean_T *y;
  emxArray_cell_wrap_6 *mat_tmp;
  emxArray_char_T *b_fid;
  emxArray_char_T *line;
  emxArray_char_T *r;
  emxArray_char_T *tok;
  emxArray_int32_T *b_n;
  emxArray_int32_T *ind2;
  emxArray_int32_T *r1;
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
  int exitg1;
  int fid;
  int i;
  int i1;
  int j;
  int k;
  int nz;
  int vlen;
  int wherefrom;
  int *ind2_data;
  int *r2;
  char a[5];
  const char *filename_data;
  signed char fileid;
  char *line_data;
  bool b_bool;
  bool *y_data;
  m_data = m->data;
  filename_data = filename->data;
  /*  filename is the meme file that contains the motifs; */
  /*  n is the nth motif in the file */
  /* 'getmotif:4' mat = cell(length(m),1); */
  vlen = m->size[0];
  i = mat->size[0];
  mat->size[0] = m->size[0];
  emxEnsureCapacity_cell_wrap_6(mat, i);
  mat_data = mat->data;
  for (i = 0; i < vlen; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 0;
  }
  emxInit_real_T(&n, 1);
  /* 'getmotif:5' mat = coder.nullcopy(mat); */
  /* 'getmotif:6' names = cell(length(m),1); */
  /* 'getmotif:7' names = coder.nullcopy(names); */
  /* 'getmotif:8' [n ind] = sort(m); */
  i = n->size[0];
  n->size[0] = m->size[0];
  emxEnsureCapacity_real_T(n, i);
  n_data = n->data;
  vlen = m->size[0];
  for (i = 0; i < vlen; i++) {
    n_data[i] = m_data[i];
  }
  emxInit_int32_T(&b_n, 1);
  sort(n, b_n);
  n_data = n->data;
  /* 'getmotif:9' fid = fopen(filename); */
  b_bool = false;
  emxFree_int32_T(&b_n);
  if (filename->size[1] == 3) {
    vlen = 0;
    do {
      exitg1 = 0;
      if (vlen < 3) {
        if (filename_data[vlen] != b_cv[vlen]) {
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
  /* 'getmotif:10' line = ''; */
  /* 'getmotif:11' if fid < 0 */
  if (fid < 0) {
    /* 'getmotif:12' fprintf("ERROR: Cannot open motif database.\n"); */
    printf("ERROR: Cannot open motif database.\n");
    fflush(stdout);
  } else {
    /* 'getmotif:13' else */
    /* 'getmotif:14' i=0; */
    b_i = 0.0;
    /* 'getmotif:15' for j = 1:length(n) */
    i = n->size[0];
    emxInit_char_T(&line, 2);
    emxInit_real_T(&zero_idx, 2);
    emxInit_real_T(&tmp, 2);
    emxInit_char_T(&tok, 2);
    emxInit_boolean_T(&y, 2);
    emxInit_char_T(&r, 2);
    emxInit_char_T(&b_fid, 2);
    for (j = 0; j < i; j++) {
      /* 'getmotif:16' while i~=n(j) */
      while (b_i != n_data[j]) {
        /* 'getmotif:17' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
        /* 'getmotif:18' if length(line) >= 5 */
        if (line->size[1] >= 5) {
          /* 'getmotif:19' if strcmp(line(1:5), 'MOTIF') */
          for (i1 = 0; i1 < 5; i1++) {
            a[i1] = line_data[i1];
          }
          vlen = memcmp(&a[0], &b[0], 5);
          if (vlen == 0) {
            /* 'getmotif:20' i = i+1; */
            b_i++;
          }
        }
      }
      /* 'getmotif:24' line_clean = strip(line); */
      /* 'getmotif:25' zero_idx = strfind(line_clean, ' '); */
      /* 'getmotif:26' target_zero = zero_idx(end); */
      /* 'getmotif:27' names{j} = line_clean(target_zero+1:end); */
      /*  a = strsplit(line, ' '); */
      /*  names{i} = a{end}; */
      /* 'getmotif:30' line = fgetl(fid); */
      fgetl(fid, b_fid);
      /* 'getmotif:32' loc = ftell(fid); */
      getfilestar(fid, &filestar, &b_bool);
      if ((fid == 0) || (fid == 1) || (fid == 2)) {
        filestar = NULL;
      }
      if (filestar == NULL) {
        loc = -1.0;
      } else {
        position_t = ftell(filestar);
        loc = (double)position_t;
      }
      /* 'getmotif:33' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;
      /* 'getmotif:35' row_count = 0; */
      row_count = 0.0;
      /* 'getmotif:36' col_count = 0; */
      col_count = 0.0;
      /* 'getmotif:37' while ~isempty(line) */
      while (line->size[1] != 0) {
        /* 'getmotif:38' row_count = row_count + 1; */
        row_count++;
        /* 'getmotif:39' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:40' col_count         = sum(diff(number_whitespace)==1) +
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
        diff(y, zero_idx);
        zero_idx_data = zero_idx->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = zero_idx->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        vlen = zero_idx->size[1];
        for (i1 = 0; i1 < vlen; i1++) {
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
        /* 'getmotif:41' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
      }
      /* 'getmotif:43' fseek(fid, loc, "bof"); */
      wherefrom = SEEK_SET;
      if ((!rtIsInf(loc)) && (!rtIsNaN(loc)) && (floor(loc) == loc)) {
        getfilestar(fid, &filestar, &b_bool);
        if ((fid == 0) || (fid == 1) || (fid == 2)) {
          filestar = NULL;
        }
        if (!(filestar == NULL)) {
          fseek(filestar, (long int)loc, wherefrom);
        }
      }
      /* 'getmotif:44' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;
      /* 'getmotif:46' tmp = zeros(row_count, col_count); */
      i1 = tmp->size[0] * tmp->size[1];
      tmp->size[0] = (int)row_count;
      tmp->size[1] = (int)col_count;
      emxEnsureCapacity_real_T(tmp, i1);
      tmp_data = tmp->data;
      vlen = (int)row_count * (int)col_count;
      for (i1 = 0; i1 < vlen; i1++) {
        tmp_data[i1] = 0.0;
      }
      /* 'getmotif:47' curr_row = 1; */
      loc = 1.0;
      /* 'getmotif:48' while ~isempty(line) */
      while (line->size[1] != 0) {
        /*  Alternative to str2num */
        /*  mat{j} = [mat{j}; str2num(line)]; */
        /* 'getmotif:51' number_whitespace = isstrprop(line, 'wspace'); */
        /* 'getmotif:52' number_element    = sum(diff(number_whitespace)==1) +
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
        diff(y, zero_idx);
        zero_idx_data = zero_idx->data;
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = zero_idx->size[1];
        emxEnsureCapacity_boolean_T(y, i1);
        y_data = y->data;
        vlen = zero_idx->size[1];
        for (i1 = 0; i1 < vlen; i1++) {
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
        /* 'getmotif:53' remain  = line; */
        /* 'getmotif:54' for idx=1:number_element */
        for (vlen = 0; vlen <= nz; vlen++) {
          /* 'getmotif:55' [tok, remain] = strtok(strip(remain), [char(32),
           * char(9)]); */
          strip(line, r);
          d_strtok(r, tok, line);
          /* 'getmotif:56' tmp(curr_row, idx)  = real(str2double(tok)); */
          dc = str2double(tok);
          tmp_data[((int)loc + tmp->size[0] * vlen) - 1] = dc.re;
        }
        /* 'getmotif:58' line = fgetl(fid); */
        fgetl(fid, line);
        line_data = line->data;
        /* 'getmotif:59' curr_row = curr_row + 1; */
        loc++;
      }
      /* 'getmotif:61' mat{j} = tmp; */
      i1 = mat_data[j].f1->size[0] * mat_data[j].f1->size[1];
      mat_data[j].f1->size[0] = tmp->size[0];
      mat_data[j].f1->size[1] = tmp->size[1];
      emxEnsureCapacity_real_T(mat_data[j].f1, i1);
      vlen = tmp->size[0] * tmp->size[1];
      for (i1 = 0; i1 < vlen; i1++) {
        mat_data[j].f1->data[i1] = tmp_data[i1];
      }
    }
    emxFree_char_T(&b_fid);
    emxFree_char_T(&r);
    emxFree_boolean_T(&y);
    emxFree_char_T(&tok);
    emxFree_real_T(&tmp);
    emxFree_real_T(&zero_idx);
    emxFree_char_T(&line);
  }
  /* 'getmotif:64' fclose(fid); */
  cfclose(fid);
  /* 'getmotif:65' if length(n)==1 */
  if (n->size[0] != 1) {
    emxInit_int32_T(&ind2, 1);
    /* 'getmotif:68' else */
    /* 'getmotif:69' mylength = length(n); */
    /* 'getmotif:70' ind2 = zeros(mylength,1); */
    i = ind2->size[0];
    ind2->size[0] = n->size[0];
    emxEnsureCapacity_int32_T(ind2, i);
    ind2_data = ind2->data;
    vlen = n->size[0];
    for (i = 0; i < vlen; i++) {
      ind2_data[i] = 0;
    }
    /* 'getmotif:71' for i = 1:mylength */
    i = n->size[0];
    for (nz = 0; nz < i; nz++) {
      /* 'getmotif:72' ind2(i,1) = i; */
      ind2_data[nz] = nz + 1;
    }
    /* 'getmotif:74' [a, b] = size(ind2); */
    /* 'getmotif:75' for i = 1:length(n) */
    i = n->size[0];
    emxInit_boolean_T(&c_n, 1);
    emxInit_int32_T(&r1, 1);
    for (nz = 0; nz < i; nz++) {
      /* 'getmotif:76' ind2(i,1) = find(n==m(i)); */
      vlen = n->size[0];
      i1 = c_n->size[0];
      c_n->size[0] = n->size[0];
      emxEnsureCapacity_boolean_T(c_n, i1);
      y_data = c_n->data;
      for (i1 = 0; i1 < vlen; i1++) {
        y_data[i1] = (n_data[i1] == m_data[nz]);
      }
      eml_find(c_n, r1);
      r2 = r1->data;
      ind2_data[nz] = r2[0];
    }
    emxFree_int32_T(&r1);
    emxFree_boolean_T(&c_n);
    emxInit_cell_wrap_6(&mat_tmp, 1);
    /* 'getmotif:78' mat_tmp = cell(mylength,1); */
    /* 'getmotif:79' for i = 1:mylength */
    i = n->size[0];
    i1 = mat_tmp->size[0];
    mat_tmp->size[0] = n->size[0];
    emxEnsureCapacity_cell_wrap_6(mat_tmp, i1);
    mat_tmp_data = mat_tmp->data;
    for (nz = 0; nz < i; nz++) {
      /* 'getmotif:80' mat_tmp{i} = mat{ind2(i)}; */
      i1 = mat_tmp_data[nz].f1->size[0] * mat_tmp_data[nz].f1->size[1];
      mat_tmp_data[nz].f1->size[0] = mat_data[ind2_data[nz] - 1].f1->size[0];
      mat_tmp_data[nz].f1->size[1] = mat_data[ind2_data[nz] - 1].f1->size[1];
      emxEnsureCapacity_real_T(mat_tmp_data[nz].f1, i1);
      vlen = mat_data[ind2_data[nz] - 1].f1->size[0] *
             mat_data[ind2_data[nz] - 1].f1->size[1];
      for (i1 = 0; i1 < vlen; i1++) {
        mat_tmp_data[nz].f1->data[i1] =
            mat_data[ind2_data[nz] - 1].f1->data[i1];
      }
    }
    emxFree_int32_T(&ind2);
    /* 'getmotif:82' [a, b] = size(mat_tmp{1}); */
    /* 'getmotif:83' for i = 1:length(n) */
    i = n->size[0];
    for (nz = 0; nz < i; nz++) {
      /* 'getmotif:84' mat{i} = mat_tmp{i}; */
      i1 = mat_data[nz].f1->size[0] * mat_data[nz].f1->size[1];
      mat_data[nz].f1->size[0] = mat_tmp_data[nz].f1->size[0];
      mat_data[nz].f1->size[1] = mat_tmp_data[nz].f1->size[1];
      emxEnsureCapacity_real_T(mat_data[nz].f1, i1);
      vlen = mat_tmp_data[nz].f1->size[0] * mat_tmp_data[nz].f1->size[1];
      for (i1 = 0; i1 < vlen; i1++) {
        mat_data[nz].f1->data[i1] = mat_tmp_data[nz].f1->data[i1];
      }
    }
    emxFree_cell_wrap_6(&mat_tmp);
  } else {
    /*  Alternative to cell2mat */
    /*  mat = cell2mat(mat); */
  }
  emxFree_real_T(&n);
}

/* End of code generation (getmotif.c) */
