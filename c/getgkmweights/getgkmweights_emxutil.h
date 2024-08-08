/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmweights_emxutil.h
 *
 * Code generation for function 'getgkmweights_emxutil'
 *
 */

#ifndef GETGKMWEIGHTS_EMXUTIL_H
#define GETGKMWEIGHTS_EMXUTIL_H

/* Include files */
#include "getgkmweights_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void emxCopyStruct_cell_wrap_6(cell_wrap_6 *dst, const cell_wrap_6 *src);

extern void emxCopy_char_T(emxArray_char_T **dst, emxArray_char_T *const *src);

extern void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray,
                                        int oldNumel);

extern void emxEnsureCapacity_cell_wrap_10(emxArray_cell_wrap_10 *emxArray,
                                           int oldNumel);

extern void emxEnsureCapacity_cell_wrap_7(emxArray_cell_wrap_7 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_uint32_T(emxArray_uint32_T *emxArray,
                                       int oldNumel);

extern void emxExpand_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxFreeMatrix_cell_wrap_1(cell_wrap_1 pMatrix[8]);

extern void emxFreeMatrix_cell_wrap_6(cell_wrap_6 *pMatrix);

extern void emxFreeMatrix_cell_wrap_61(cell_wrap_6 pMatrix[3]);

extern void emxFreeStruct_cell_wrap_1(cell_wrap_1 *pStruct);

extern void emxFreeStruct_cell_wrap_6(cell_wrap_6 *pStruct);

extern void emxFreeStruct_cell_wrap_7(cell_wrap_7 *pStruct);

extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);

extern void emxFree_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray);

extern void emxFree_cell_wrap_7(emxArray_cell_wrap_7 **pEmxArray);

extern void emxFree_char_T(emxArray_char_T **pEmxArray);

extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);

extern void emxFree_real_T(emxArray_real_T **pEmxArray);

extern void emxFree_uint32_T(emxArray_uint32_T **pEmxArray);

extern void emxInitMatrix_cell_wrap_1(cell_wrap_1 pMatrix[8]);

extern void emxInitMatrix_cell_wrap_6(cell_wrap_6 *pMatrix);

extern void emxInitMatrix_cell_wrap_61(cell_wrap_6 pMatrix[3]);

extern void emxInitStruct_cell_wrap_1(cell_wrap_1 *pStruct);

extern void emxInitStruct_cell_wrap_6(cell_wrap_6 *pStruct);

extern void emxInitStruct_cell_wrap_7(cell_wrap_7 *pStruct);

extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray,
                              int numDimensions);

extern void emxInit_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray);

extern void emxInit_cell_wrap_7(emxArray_cell_wrap_7 **pEmxArray);

extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);

extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);

extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

extern void emxInit_uint32_T(emxArray_uint32_T **pEmxArray);

extern void emxTrim_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int fromIndex,
                                int toIndex);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (getgkmweights_emxutil.h) */
