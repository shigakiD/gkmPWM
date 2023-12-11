/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * graph.c
 *
 * Code generation for function 'graph'
 *
 */

/* Include files */
#include "graph.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void graph_graph(const emxArray_boolean_T *varargin_1,
                 emxArray_int32_T *G_Underlying_Ir,
                 emxArray_int32_T *G_Underlying_Jc)
{
  emxArray_int32_T *b_i;
  emxArray_int32_T *b_r;
  emxArray_int32_T *irNew;
  emxArray_int32_T *j;
  emxArray_int32_T *listOfSelfLoops;
  emxArray_int32_T *r2;
  int i;
  int i1;
  int i2;
  int idx;
  int ii;
  int jj;
  int nx;
  int *i_data;
  int *irNew_data;
  int *j_data;
  int *listOfSelfLoops_data;
  int *r1;
  int *r3;
  const bool *varargin_1_data;
  bool exitg1;
  bool guard1 = false;
  varargin_1_data = varargin_1->data;
  nx = varargin_1->size[0] * varargin_1->size[1];
  emxInit_int32_T(&j, 1);
  j_data = j->data;
  emxInit_int32_T(&irNew, 1);
  irNew_data = irNew->data;
  if (nx == 0) {
    irNew->size[0] = 0;
    j->size[0] = 0;
  } else {
    idx = 0;
    i = irNew->size[0];
    irNew->size[0] = nx;
    emxEnsureCapacity_int32_T(irNew, i);
    irNew_data = irNew->data;
    i = j->size[0];
    j->size[0] = nx;
    emxEnsureCapacity_int32_T(j, i);
    j_data = j->data;
    ii = 1;
    jj = 1;
    exitg1 = false;
    while ((!exitg1) && (jj <= varargin_1->size[1])) {
      guard1 = false;
      if (varargin_1_data[(ii + varargin_1->size[0] * (jj - 1)) - 1]) {
        idx++;
        irNew_data[idx - 1] = ii;
        j_data[idx - 1] = jj;
        if (idx >= nx) {
          exitg1 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
      if (guard1) {
        ii++;
        if (ii > varargin_1->size[0]) {
          ii = 1;
          jj++;
        }
      }
    }
    if (nx == 1) {
      if (idx == 0) {
        irNew->size[0] = 0;
        j->size[0] = 0;
      }
    } else {
      i = irNew->size[0];
      if (1 > idx) {
        irNew->size[0] = 0;
      } else {
        irNew->size[0] = idx;
      }
      emxEnsureCapacity_int32_T(irNew, i);
      irNew_data = irNew->data;
      i = j->size[0];
      if (1 > idx) {
        j->size[0] = 0;
      } else {
        j->size[0] = idx;
      }
      emxEnsureCapacity_int32_T(j, i);
      j_data = j->data;
    }
  }
  emxInit_int32_T(&b_i, 1);
  i = b_i->size[0];
  b_i->size[0] = irNew->size[0];
  emxEnsureCapacity_int32_T(b_i, i);
  i_data = b_i->data;
  jj = irNew->size[0];
  for (i = 0; i < jj; i++) {
    i_data[i] = irNew_data[i];
  }
  i = G_Underlying_Jc->size[0];
  G_Underlying_Jc->size[0] = varargin_1->size[0] + 1;
  emxEnsureCapacity_int32_T(G_Underlying_Jc, i);
  irNew_data = G_Underlying_Jc->data;
  jj = varargin_1->size[0];
  for (i = 0; i <= jj; i++) {
    irNew_data[i] = 0;
  }
  emxInit_int32_T(&b_r, 2);
  i = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = b_i->size[0];
  emxEnsureCapacity_int32_T(b_r, i);
  r1 = b_r->data;
  jj = b_i->size[0];
  for (i = 0; i < jj; i++) {
    r1[i] = i_data[i];
  }
  jj = 0;
  i = varargin_1->size[0];
  for (ii = 0; ii < i; ii++) {
    while ((jj < b_i->size[0]) && (j_data[jj] < (double)ii + 2.0)) {
      jj++;
    }
    irNew_data[ii + 1] = jj;
  }
  emxFree_int32_T(&j);
  i = G_Underlying_Ir->size[0] * G_Underlying_Ir->size[1];
  G_Underlying_Ir->size[0] = 1;
  G_Underlying_Ir->size[1] = b_i->size[0];
  emxEnsureCapacity_int32_T(G_Underlying_Ir, i);
  j_data = G_Underlying_Ir->data;
  jj = b_i->size[0];
  for (i = 0; i < jj; i++) {
    j_data[i] = i_data[i];
  }
  emxInit_int32_T(&r2, 2);
  emxInit_int32_T(&listOfSelfLoops, 2);
  listOfSelfLoops_data = listOfSelfLoops->data;
  listOfSelfLoops->size[0] = 1;
  listOfSelfLoops->size[1] = 0;
  i = r2->size[0] * r2->size[1];
  r2->size[0] = 1;
  r2->size[1] = b_i->size[0];
  emxEnsureCapacity_int32_T(r2, i);
  r3 = r2->data;
  jj = b_i->size[0];
  for (i = 0; i < jj; i++) {
    r3[i] = i_data[i];
  }
  i = varargin_1->size[0];
  for (nx = 0; nx < i; nx++) {
    jj = listOfSelfLoops->size[1];
    ii = irNew_data[nx];
    exitg1 = false;
    while ((!exitg1) && (ii + 1 <= irNew_data[nx + 1])) {
      if (nx + 1 == j_data[ii]) {
        i1 = listOfSelfLoops->size[1];
        i2 = listOfSelfLoops->size[0] * listOfSelfLoops->size[1];
        listOfSelfLoops->size[1]++;
        emxEnsureCapacity_int32_T(listOfSelfLoops, i2);
        listOfSelfLoops_data = listOfSelfLoops->data;
        listOfSelfLoops_data[i1] = ii + 1;
        ii++;
      } else if (j_data[ii] > nx + 1) {
        exitg1 = true;
      } else {
        ii++;
      }
    }
    irNew_data[nx] += jj;
  }
  irNew_data[G_Underlying_Jc->size[0] - 1] += listOfSelfLoops->size[1];
  nx = 0;
  idx = 1;
  if (listOfSelfLoops->size[1] != 0) {
    if (G_Underlying_Ir->size[1] == 0) {
      jj = 0;
    } else {
      jj = G_Underlying_Ir->size[1];
    }
    i = irNew->size[0];
    irNew->size[0] = jj + listOfSelfLoops->size[1];
    emxEnsureCapacity_int32_T(irNew, i);
    irNew_data = irNew->data;
    jj += listOfSelfLoops->size[1];
    for (i = 0; i < jj; i++) {
      irNew_data[i] = 0;
    }
    i = listOfSelfLoops->size[1];
    for (ii = 0; ii < i; ii++) {
      i1 = listOfSelfLoops_data[ii];
      if (idx > i1) {
        i2 = 0;
        i1 = 0;
      } else {
        i2 = idx - 1;
      }
      jj = i1 - i2;
      for (i1 = 0; i1 < jj; i1++) {
        irNew_data[i1 + nx] = r1[i2 + i1];
      }
      i1 = listOfSelfLoops_data[ii];
      nx = (nx + i1) - idx;
      idx = i1;
      irNew_data[nx] = r1[i1 - 1];
      nx++;
    }
    if (idx > b_i->size[0]) {
      i = 0;
      i1 = 0;
    } else {
      i = idx - 1;
      i1 = G_Underlying_Ir->size[1];
    }
    if (nx + 1 > irNew->size[0]) {
      nx = 0;
    }
    jj = i1 - i;
    for (i1 = 0; i1 < jj; i1++) {
      irNew_data[nx + i1] = r1[i + i1];
    }
    i = r2->size[0] * r2->size[1];
    r2->size[0] = irNew->size[0];
    r2->size[1] = 1;
    emxEnsureCapacity_int32_T(r2, i);
    r3 = r2->data;
    jj = irNew->size[0];
    for (i = 0; i < jj; i++) {
      r3[i] = irNew_data[i];
    }
  }
  emxFree_int32_T(&irNew);
  emxFree_int32_T(&listOfSelfLoops);
  emxFree_int32_T(&b_r);
  emxFree_int32_T(&b_i);
  i = G_Underlying_Ir->size[0] * G_Underlying_Ir->size[1];
  G_Underlying_Ir->size[0] = r2->size[0];
  G_Underlying_Ir->size[1] = r2->size[1];
  emxEnsureCapacity_int32_T(G_Underlying_Ir, i);
  j_data = G_Underlying_Ir->data;
  jj = r2->size[0] * r2->size[1];
  for (i = 0; i < jj; i++) {
    j_data[i] = r3[i];
  }
  emxFree_int32_T(&r2);
}

/* End of code generation (graph.c) */
