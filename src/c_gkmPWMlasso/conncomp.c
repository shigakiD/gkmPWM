/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * conncomp.c
 *
 * Code generation for function 'conncomp'
 *
 */

/* Include files */
#include "conncomp.h"
#include "gkmPWMlasso4_data.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include "list.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void graph_conncomp(const emxArray_int32_T *G_Underlying_Ir,
                    const emxArray_int32_T *G_Underlying_Jc,
                    emxArray_real_T *bins)
{
  coder_internal_list nodeList;
  double nextbin;
  double *bins_data;
  const int *G_Underlying_Ir_data;
  const int *G_Underlying_Jc_data;
  int i;
  int i1;
  int k;
  int l;
  int s;
  int start;
  G_Underlying_Jc_data = G_Underlying_Jc->data;
  G_Underlying_Ir_data = G_Underlying_Ir->data;
  i = bins->size[0] * bins->size[1];
  bins->size[0] = 1;
  bins->size[1] = G_Underlying_Jc->size[0] - 1;
  emxEnsureCapacity_real_T(bins, i);
  bins_data = bins->data;
  k = G_Underlying_Jc->size[0] - 1;
  for (i = 0; i < k; i++) {
    bins_data[i] = 0.0;
  }
  nextbin = 0.0;
  i = G_Underlying_Jc->size[0];
  c_emxInitStruct_coder_internal_(&nodeList);
  for (start = 0; start <= i - 2; start++) {
    if (bins_data[start] == 0.0) {
      nextbin++;
      k = nodeList.valuePool->size[0] * nodeList.valuePool->size[1];
      nodeList.valuePool->size[0] = 1;
      nodeList.valuePool->size[1] = 1;
      emxEnsureCapacity_int32_T(nodeList.valuePool, k);
      nodeList.valuePool->data[0] = start + 1;
      k = nodeList.nodePool->size[0] * nodeList.nodePool->size[1];
      nodeList.nodePool->size[0] = 1;
      nodeList.nodePool->size[1] = 1;
      emxEnsureCapacity_struct_T(nodeList.nodePool, k);
      nodeList.nodePool->data[0] = r;
      nodeList.len = 0;
      nodeList.frontAddr = 0;
      nodeList.backAddr = 0;
      nodeList.nodePool->data[0].addr = 1;
      nodeList.nodePool->data[0].next = 2;
      nodeList.nodePool->data[0].prev = 0;
      nodeList.nodePool->data[0].next = 0;
      nodeList.unusedAddr = 1;
      list_pushBack(&nodeList, start + 1);
      while (nodeList.len > 0) {
        if (nodeList.frontAddr == 0) {
          s = 0;
        } else {
          k = nodeList.frontAddr;
          s = nodeList.valuePool->data[nodeList.frontAddr - 1];
          nodeList.frontAddr =
              nodeList.nodePool->data[nodeList.frontAddr - 1].next;
          if (nodeList.frontAddr != 0) {
            nodeList.nodePool->data[nodeList.frontAddr - 1].prev = 0;
          } else {
            nodeList.backAddr = 0;
          }
          nodeList.len--;
          if (nodeList.unusedAddr != 0) {
            nodeList.nodePool->data[nodeList.unusedAddr - 1].prev = k;
          }
          nodeList.nodePool->data[k - 1].next = nodeList.unusedAddr;
          nodeList.nodePool->data[k - 1].prev = 0;
          nodeList.unusedAddr = k;
        }
        k = G_Underlying_Jc_data[s - 1] + 1;
        i1 = G_Underlying_Jc_data[s];
        for (l = k; l <= i1; l++) {
          if (bins_data[G_Underlying_Ir_data[l - 1] - 1] == 0.0) {
            bins_data[G_Underlying_Ir_data[l - 1] - 1] = -1.0;
            list_pushBack(&nodeList, G_Underlying_Ir_data[l - 1]);
          }
        }
        bins_data[s - 1] = nextbin;
      }
    }
  }
  c_emxFreeStruct_coder_internal_(&nodeList);
}

/* End of code generation (conncomp.c) */
