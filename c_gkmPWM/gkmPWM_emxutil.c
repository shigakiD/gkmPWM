/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM_emxutil.c
 *
 * Code generation for function 'gkmPWM_emxutil'
 *
 */

/* Include files */
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "lapacke.h"
#include <stdlib.h>
#include <string.h>

/* Function Definitions */
void emxCopyStruct_cell_wrap_0(cell_wrap_0 *dst, const cell_wrap_0 *src)
{
  emxCopy_char_T(&dst->f1, &src->f1);
}

void emxCopyStruct_cell_wrap_2(cell_wrap_2 *dst, const cell_wrap_2 *src)
{
  emxCopy_real_T(&dst->f1, &src->f1);
}

void emxCopy_char_T(emxArray_char_T **dst, emxArray_char_T *const *src)
{
  int i;
  int numElDst;
  int numElSrc;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }
  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }
  emxEnsureCapacity_char_T(*dst, numElDst);
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T *const *src)
{
  int i;
  int numElDst;
  int numElSrc;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }
  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }
  emxEnsureCapacity_real_T(*dst, numElDst);
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(bool));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(bool) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (bool *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(cell_wrap_0));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_0) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_0 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_cell_wrap_0(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_cell_wrap_0(emxArray, oldNumel, newNumel);
  }
}

void emxEnsureCapacity_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(cell_wrap_1));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_1) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_1 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_cell_wrap_1(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_cell_wrap_1(emxArray, oldNumel, newNumel);
  }
}

void emxEnsureCapacity_cell_wrap_10(emxArray_cell_wrap_10 *emxArray,
                                    int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(cell_wrap_10));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_10) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_10 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_cell_wrap_10(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_cell_wrap_10(emxArray, oldNumel, newNumel);
  }
}

void emxEnsureCapacity_cell_wrap_13(emxArray_cell_wrap_13 *emxArray,
                                    int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(cell_wrap_13));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_13) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_13 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_cell_wrap_14(emxArray_cell_wrap_14 *emxArray,
                                    int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(cell_wrap_14));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_14) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_14 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_cell_wrap_14(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_cell_wrap_14(emxArray, oldNumel, newNumel);
  }
}

void emxEnsureCapacity_cell_wrap_2(emxArray_cell_wrap_2 *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(cell_wrap_2));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_2) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_2 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_cell_wrap_2(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_cell_wrap_2(emxArray, oldNumel, newNumel);
  }
}

void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(char));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(char) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (char *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_creal_T(emxArray_creal_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(creal_T));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(creal_T) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (creal_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(int));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(int) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (int *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(signed char));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(signed char) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (signed char *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_lapack_int(emxArray_lapack_int *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(lapack_int));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(lapack_int) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (lapack_int *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(double));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(double) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (double *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxEnsureCapacity_rtString(emxArray_rtString *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(rtString));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(rtString) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (rtString *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_rtString(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_rtString(emxArray, oldNumel, newNumel);
  }
}

void emxEnsureCapacity_uint32_T(emxArray_uint32_T *emxArray, int oldNumel)
{
  int i;
  int newNumel;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }
  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }
    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }
    newData = calloc((unsigned int)i, sizeof(unsigned int));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(unsigned int) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (unsigned int *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

void emxExpand_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int fromIndex,
                           int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_cell_wrap_0(&emxArray->data[i]);
  }
}

void emxExpand_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int fromIndex,
                           int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_cell_wrap_1(&emxArray->data[i]);
  }
}

void emxExpand_cell_wrap_10(emxArray_cell_wrap_10 *emxArray, int fromIndex,
                            int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_cell_wrap_10(&emxArray->data[i]);
  }
}

void emxExpand_cell_wrap_14(emxArray_cell_wrap_14 *emxArray, int fromIndex,
                            int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_cell_wrap_14(&emxArray->data[i]);
  }
}

void emxExpand_cell_wrap_2(emxArray_cell_wrap_2 *emxArray, int fromIndex,
                           int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_cell_wrap_2(&emxArray->data[i]);
  }
}

void emxExpand_rtString(emxArray_rtString *emxArray, int fromIndex, int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_rtString(&emxArray->data[i]);
  }
}

void emxFreeMatrix_cell_wrap_0(cell_wrap_0 *pMatrix)
{
  emxFreeStruct_cell_wrap_0(pMatrix);
}

void emxFreeMatrix_cell_wrap_3(cell_wrap_3 pMatrix[1968])
{
  int i;
  for (i = 0; i < 1968; i++) {
    emxFreeStruct_cell_wrap_3(&pMatrix[i]);
  }
}

void emxFreeMatrix_cell_wrap_31(cell_wrap_3 pMatrix[1969])
{
  int i;
  for (i = 0; i < 1969; i++) {
    emxFreeStruct_cell_wrap_3(&pMatrix[i]);
  }
}

void emxFreeStruct_cell_wrap_0(cell_wrap_0 *pStruct)
{
  emxFree_char_T(&pStruct->f1);
}

void emxFreeStruct_cell_wrap_1(cell_wrap_1 *pStruct)
{
  emxFree_real_T(&pStruct->f1);
}

void emxFreeStruct_cell_wrap_10(cell_wrap_10 *pStruct)
{
  emxFreeMatrix_cell_wrap_0(&pStruct->f1[0]);
}

void emxFreeStruct_cell_wrap_14(cell_wrap_14 *pStruct)
{
  emxFree_real_T(&pStruct->f1);
}

void emxFreeStruct_cell_wrap_2(cell_wrap_2 *pStruct)
{
  emxFree_real_T(&pStruct->f1);
}

void emxFreeStruct_cell_wrap_3(cell_wrap_3 *pStruct)
{
  emxFree_real_T(&pStruct->f1);
}

void emxFreeStruct_rtString(rtString *pStruct)
{
  emxFree_char_T(&pStruct->Value);
}

void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (bool *)NULL) && (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

void emxFree_cell_wrap_0(emxArray_cell_wrap_0 **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_cell_wrap_0 *)NULL) {
    if ((*pEmxArray)->data != (cell_wrap_0 *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_cell_wrap_0(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_0 *)NULL;
  }
}

void emxFree_cell_wrap_1(emxArray_cell_wrap_1 **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_cell_wrap_1 *)NULL) {
    if ((*pEmxArray)->data != (cell_wrap_1 *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_cell_wrap_1(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_1 *)NULL;
  }
}

void emxFree_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_cell_wrap_10 *)NULL) {
    if ((*pEmxArray)->data != (cell_wrap_10 *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_cell_wrap_10(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_10 *)NULL;
  }
}

void emxFree_cell_wrap_13(emxArray_cell_wrap_13 **pEmxArray)
{
  if (*pEmxArray != (emxArray_cell_wrap_13 *)NULL) {
    if (((*pEmxArray)->data != (cell_wrap_13 *)NULL) &&
        (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_13 *)NULL;
  }
}

void emxFree_cell_wrap_14(emxArray_cell_wrap_14 **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_cell_wrap_14 *)NULL) {
    if ((*pEmxArray)->data != (cell_wrap_14 *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_cell_wrap_14(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_14 *)NULL;
  }
}

void emxFree_cell_wrap_2(emxArray_cell_wrap_2 **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_cell_wrap_2 *)NULL) {
    if ((*pEmxArray)->data != (cell_wrap_2 *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_cell_wrap_2(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_2 *)NULL;
  }
}

void emxFree_char_T(emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if (((*pEmxArray)->data != (char *)NULL) && (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

void emxFree_creal_T(emxArray_creal_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_creal_T *)NULL) {
    if (((*pEmxArray)->data != (creal_T *)NULL) && (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_creal_T *)NULL;
  }
}

void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

void emxFree_int8_T(emxArray_int8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int8_T *)NULL) {
    if (((*pEmxArray)->data != (signed char *)NULL) &&
        (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_int8_T *)NULL;
  }
}

void emxFree_lapack_int(emxArray_lapack_int **pEmxArray)
{
  if (*pEmxArray != (emxArray_lapack_int *)NULL) {
    if (((*pEmxArray)->data != (lapack_int *)NULL) &&
        (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_lapack_int *)NULL;
  }
}

void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

void emxFree_rtString(emxArray_rtString **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_rtString *)NULL) {
    if ((*pEmxArray)->data != (rtString *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_rtString(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_rtString *)NULL;
  }
}

void emxFree_uint32_T(emxArray_uint32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint32_T *)NULL) {
    if (((*pEmxArray)->data != (unsigned int *)NULL) &&
        (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_uint32_T *)NULL;
  }
}

void emxInitMatrix_cell_wrap_0(cell_wrap_0 *pMatrix)
{
  emxInitStruct_cell_wrap_0(pMatrix);
}

void emxInitMatrix_cell_wrap_3(cell_wrap_3 pMatrix[1968])
{
  int i;
  for (i = 0; i < 1968; i++) {
    emxInitStruct_cell_wrap_3(&pMatrix[i]);
  }
}

void emxInitMatrix_cell_wrap_31(cell_wrap_3 pMatrix[1969])
{
  int i;
  for (i = 0; i < 1969; i++) {
    emxInitStruct_cell_wrap_3(&pMatrix[i]);
  }
}

void emxInitStruct_cell_wrap_0(cell_wrap_0 *pStruct)
{
  emxInit_char_T(&pStruct->f1, 2);
}

void emxInitStruct_cell_wrap_1(cell_wrap_1 *pStruct)
{
  emxInit_real_T(&pStruct->f1, 2);
}

void emxInitStruct_cell_wrap_10(cell_wrap_10 *pStruct)
{
  emxInitMatrix_cell_wrap_0(&pStruct->f1[0]);
}

void emxInitStruct_cell_wrap_14(cell_wrap_14 *pStruct)
{
  emxInit_real_T(&pStruct->f1, 1);
}

void emxInitStruct_cell_wrap_2(cell_wrap_2 *pStruct)
{
  emxInit_real_T(&pStruct->f1, 2);
}

void emxInitStruct_cell_wrap_3(cell_wrap_3 *pStruct)
{
  emxInit_real_T(&pStruct->f1, 2);
}

void emxInitStruct_rtString(rtString *pStruct)
{
  emxInit_char_T(&pStruct->Value, 2);
}

void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
  emxArray_boolean_T *emxArray;
  int i;
  *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (bool *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_cell_wrap_0(emxArray_cell_wrap_0 **pEmxArray, int numDimensions)
{
  emxArray_cell_wrap_0 *emxArray;
  int i;
  *pEmxArray = (emxArray_cell_wrap_0 *)malloc(sizeof(emxArray_cell_wrap_0));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_0 *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_cell_wrap_1(emxArray_cell_wrap_1 **pEmxArray)
{
  emxArray_cell_wrap_1 *emxArray;
  *pEmxArray = (emxArray_cell_wrap_1 *)malloc(sizeof(emxArray_cell_wrap_1));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_1 *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray)
{
  emxArray_cell_wrap_10 *emxArray;
  *pEmxArray = (emxArray_cell_wrap_10 *)malloc(sizeof(emxArray_cell_wrap_10));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_10 *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_cell_wrap_13(emxArray_cell_wrap_13 **pEmxArray)
{
  emxArray_cell_wrap_13 *emxArray;
  *pEmxArray = (emxArray_cell_wrap_13 *)malloc(sizeof(emxArray_cell_wrap_13));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_13 *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_cell_wrap_14(emxArray_cell_wrap_14 **pEmxArray)
{
  emxArray_cell_wrap_14 *emxArray;
  *pEmxArray = (emxArray_cell_wrap_14 *)malloc(sizeof(emxArray_cell_wrap_14));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_14 *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_cell_wrap_2(emxArray_cell_wrap_2 **pEmxArray)
{
  emxArray_cell_wrap_2 *emxArray;
  *pEmxArray = (emxArray_cell_wrap_2 *)malloc(sizeof(emxArray_cell_wrap_2));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_2 *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxArray_char_T *emxArray;
  int i;
  *pEmxArray = (emxArray_char_T *)malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions)
{
  emxArray_creal_T *emxArray;
  int i;
  *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
  emxArray = *pEmxArray;
  emxArray->data = (creal_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions)
{
  emxArray_int8_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int8_T *)malloc(sizeof(emxArray_int8_T));
  emxArray = *pEmxArray;
  emxArray->data = (signed char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_lapack_int(emxArray_lapack_int **pEmxArray)
{
  emxArray_lapack_int *emxArray;
  *pEmxArray = (emxArray_lapack_int *)malloc(sizeof(emxArray_lapack_int));
  emxArray = *pEmxArray;
  emxArray->data = (lapack_int *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void emxInit_rtString(emxArray_rtString **pEmxArray)
{
  emxArray_rtString *emxArray;
  *pEmxArray = (emxArray_rtString *)malloc(sizeof(emxArray_rtString));
  emxArray = *pEmxArray;
  emxArray->data = (rtString *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxInit_uint32_T(emxArray_uint32_T **pEmxArray)
{
  emxArray_uint32_T *emxArray;
  *pEmxArray = (emxArray_uint32_T *)malloc(sizeof(emxArray_uint32_T));
  emxArray = *pEmxArray;
  emxArray->data = (unsigned int *)NULL;
  emxArray->numDimensions = 1;
  emxArray->size = (int *)malloc(sizeof(int));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  emxArray->size[0] = 0;
}

void emxTrim_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int fromIndex,
                         int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_cell_wrap_0(&emxArray->data[i]);
  }
}

void emxTrim_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int fromIndex,
                         int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_cell_wrap_1(&emxArray->data[i]);
  }
}

void emxTrim_cell_wrap_10(emxArray_cell_wrap_10 *emxArray, int fromIndex,
                          int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_cell_wrap_10(&emxArray->data[i]);
  }
}

void emxTrim_cell_wrap_14(emxArray_cell_wrap_14 *emxArray, int fromIndex,
                          int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_cell_wrap_14(&emxArray->data[i]);
  }
}

void emxTrim_cell_wrap_2(emxArray_cell_wrap_2 *emxArray, int fromIndex,
                         int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_cell_wrap_2(&emxArray->data[i]);
  }
}

void emxTrim_rtString(emxArray_rtString *emxArray, int fromIndex, int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_rtString(&emxArray->data[i]);
  }
}

/* End of code generation (gkmPWM_emxutil.c) */
