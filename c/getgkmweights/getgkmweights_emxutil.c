/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmweights_emxutil.c
 *
 * Code generation for function 'getgkmweights_emxutil'
 *
 */

/* Include files */
#include "getgkmweights_emxutil.h"
#include "getgkmweights_types.h"
#include <stdlib.h>
#include <string.h>

/* Function Definitions */
void emxCopyStruct_cell_wrap_6(cell_wrap_6 *dst, const cell_wrap_6 *src)
{
  emxCopy_char_T(&dst->f1, &src->f1);
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
}

void emxEnsureCapacity_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int oldNumel)
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
    newData = calloc((unsigned int)i, sizeof(cell_wrap_7));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(cell_wrap_7) * oldNumel);
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }
    emxArray->data = (cell_wrap_7 *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
  if (oldNumel > newNumel) {
    emxTrim_cell_wrap_7(emxArray, newNumel, oldNumel);
  } else if (oldNumel < newNumel) {
    emxExpand_cell_wrap_7(emxArray, oldNumel, newNumel);
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

void emxExpand_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int fromIndex,
                           int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_cell_wrap_7(&emxArray->data[i]);
  }
}

void emxFreeMatrix_cell_wrap_1(cell_wrap_1 pMatrix[8])
{
  int i;
  for (i = 0; i < 8; i++) {
    emxFreeStruct_cell_wrap_1(&pMatrix[i]);
  }
}

void emxFreeMatrix_cell_wrap_6(cell_wrap_6 *pMatrix)
{
  emxFreeStruct_cell_wrap_6(pMatrix);
}

void emxFreeMatrix_cell_wrap_61(cell_wrap_6 pMatrix[3])
{
  int i;
  for (i = 0; i < 3; i++) {
    emxFreeStruct_cell_wrap_6(&pMatrix[i]);
  }
}

void emxFreeStruct_cell_wrap_1(cell_wrap_1 *pStruct)
{
  emxFree_real_T(&pStruct->f1);
}

void emxFreeStruct_cell_wrap_6(cell_wrap_6 *pStruct)
{
  emxFree_char_T(&pStruct->f1);
}

void emxFreeStruct_cell_wrap_7(cell_wrap_7 *pStruct)
{
  emxFreeMatrix_cell_wrap_6(&pStruct->f1[0]);
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

void emxFree_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray)
{
  if (*pEmxArray != (emxArray_cell_wrap_10 *)NULL) {
    if (((*pEmxArray)->data != (cell_wrap_10 *)NULL) &&
        (*pEmxArray)->canFreeData) {
      free((*pEmxArray)->data);
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_10 *)NULL;
  }
}

void emxFree_cell_wrap_7(emxArray_cell_wrap_7 **pEmxArray)
{
  int i;
  int numEl;
  if (*pEmxArray != (emxArray_cell_wrap_7 *)NULL) {
    if ((*pEmxArray)->data != (cell_wrap_7 *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }
      for (i = 0; i < numEl; i++) {
        emxFreeStruct_cell_wrap_7(&(*pEmxArray)->data[i]);
      }
      if ((*pEmxArray)->canFreeData) {
        free((*pEmxArray)->data);
      }
    }
    free((*pEmxArray)->size);
    free(*pEmxArray);
    *pEmxArray = (emxArray_cell_wrap_7 *)NULL;
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

void emxInitMatrix_cell_wrap_1(cell_wrap_1 pMatrix[8])
{
  int i;
  for (i = 0; i < 8; i++) {
    emxInitStruct_cell_wrap_1(&pMatrix[i]);
  }
}

void emxInitMatrix_cell_wrap_6(cell_wrap_6 *pMatrix)
{
  emxInitStruct_cell_wrap_6(pMatrix);
}

void emxInitMatrix_cell_wrap_61(cell_wrap_6 pMatrix[3])
{
  int i;
  for (i = 0; i < 3; i++) {
    emxInitStruct_cell_wrap_6(&pMatrix[i]);
  }
}

void emxInitStruct_cell_wrap_1(cell_wrap_1 *pStruct)
{
  emxInit_real_T(&pStruct->f1, 2);
}

void emxInitStruct_cell_wrap_6(cell_wrap_6 *pStruct)
{
  emxInit_char_T(&pStruct->f1, 2);
}

void emxInitStruct_cell_wrap_7(cell_wrap_7 *pStruct)
{
  emxInitMatrix_cell_wrap_6(&pStruct->f1[0]);
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

void emxInit_cell_wrap_7(emxArray_cell_wrap_7 **pEmxArray)
{
  emxArray_cell_wrap_7 *emxArray;
  *pEmxArray = (emxArray_cell_wrap_7 *)malloc(sizeof(emxArray_cell_wrap_7));
  emxArray = *pEmxArray;
  emxArray->data = (cell_wrap_7 *)NULL;
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

void emxTrim_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int fromIndex,
                         int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_cell_wrap_7(&emxArray->data[i]);
  }
}

/* End of code generation (getgkmweights_emxutil.c) */
