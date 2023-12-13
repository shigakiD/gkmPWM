/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * list.c
 *
 * Code generation for function 'list'
 *
 */

/* Include files */
#include "list.h"
#include "gkmPWMlasso_data.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void list_pushBack(coder_internal_list *obj, int value)
{
  emxArray_int32_T *b_r;
  emxArray_struct_T *r2;
  struct_T *r3;
  int cap;
  int i;
  int n;
  int *r1;
  if (obj->unusedAddr == 0) {
    emxInit_int32_T(&b_r, 1);
    if (obj->nodePool->size[0] == 0) {
      n = -1;
    } else if (obj->nodePool->size[0] > 1) {
      n = obj->nodePool->size[0] - 1;
    } else {
      n = 0;
    }
    i = b_r->size[0];
    b_r->size[0] = obj->valuePool->size[0] + 8;
    emxEnsureCapacity_int32_T(b_r, i);
    r1 = b_r->data;
    cap = obj->valuePool->size[0];
    for (i = 0; i < cap; i++) {
      r1[i] = obj->valuePool->data[i];
    }
    for (i = 0; i < 8; i++) {
      r1[i + obj->valuePool->size[0]] = 0;
    }
    i = obj->valuePool->size[0] * obj->valuePool->size[1];
    obj->valuePool->size[0] = b_r->size[0];
    obj->valuePool->size[1] = 1;
    emxEnsureCapacity_int32_T(obj->valuePool, i);
    cap = b_r->size[0];
    for (i = 0; i < cap; i++) {
      obj->valuePool->data[i] = r1[i];
    }
    emxFree_int32_T(&b_r);
    emxInit_struct_T(&r2, 1);
    i = r2->size[0];
    r2->size[0] = obj->nodePool->size[0] + 8;
    emxEnsureCapacity_struct_T(r2, i);
    r3 = r2->data;
    cap = obj->nodePool->size[0];
    for (i = 0; i < cap; i++) {
      r3[i] = obj->nodePool->data[i];
    }
    for (i = 0; i < 8; i++) {
      r3[i + obj->nodePool->size[0]] = r;
    }
    i = obj->nodePool->size[0] * obj->nodePool->size[1];
    obj->nodePool->size[0] = r2->size[0];
    obj->nodePool->size[1] = 1;
    emxEnsureCapacity_struct_T(obj->nodePool, i);
    cap = r2->size[0];
    for (i = 0; i < cap; i++) {
      obj->nodePool->data[i] = r3[i];
    }
    emxFree_struct_T(&r2);
    for (cap = 0; cap < 8; cap++) {
      obj->nodePool->data[(n + cap) + 1].addr = (n + cap) + 2;
      obj->nodePool->data[(n + cap) + 1].next = (n + cap) + 3;
      obj->nodePool->data[(n + cap) + 1].prev = (n + cap) + 1;
    }
    cap = obj->nodePool->size[0];
    obj->unusedAddr = n + 2;
    obj->nodePool->data[n + 1].prev = 0;
    obj->nodePool->data[cap - 1].next = 0;
  }
  cap = obj->unusedAddr;
  if (obj->unusedAddr != 0) {
    obj->len++;
    obj->unusedAddr = obj->nodePool->data[obj->unusedAddr - 1].next;
    obj->nodePool->data[cap - 1].next = 0;
    if (obj->unusedAddr != 0) {
      obj->nodePool->data[obj->unusedAddr - 1].prev = 0;
    }
  }
  if (cap != 0) {
    if (obj->frontAddr == 0) {
      obj->frontAddr = cap;
      obj->backAddr = cap;
      obj->nodePool->data[cap - 1].next = 0;
    } else {
      obj->nodePool->data[cap - 1].prev = obj->backAddr;
      obj->nodePool->data[obj->backAddr - 1].next = cap;
      obj->backAddr = cap;
    }
    obj->valuePool->data[cap - 1] = value;
  }
}

/* End of code generation (list.c) */
