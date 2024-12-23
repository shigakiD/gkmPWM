/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF_emxutil.h
 *
 * Code generation for function 'mapTF_emxutil'
 *
 */

#ifndef MAPTF_EMXUTIL_H
#define MAPTF_EMXUTIL_H

/* Include files */
#include "mapTF_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void emxCopyStruct_cell_wrap_6(cell_wrap_6 *dst, const cell_wrap_6 *src);

extern void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T *const *src);

extern void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray,
                                        int oldNumel);

extern void emxEnsureCapacity_cell_wrap_0(emxArray_cell_wrap_0 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_1(emxArray_cell_wrap_1 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_2(emxArray_cell_wrap_2 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_3(emxArray_cell_wrap_3 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_4(emxArray_cell_wrap_4 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_5(emxArray_cell_wrap_5 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_6(emxArray_cell_wrap_6 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_7(emxArray_cell_wrap_7 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);

extern void emxExpand_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_2(emxArray_cell_wrap_2 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_3(emxArray_cell_wrap_3 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_4(emxArray_cell_wrap_4 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_5(emxArray_cell_wrap_5 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_6(emxArray_cell_wrap_6 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxFreeStruct_cell_wrap_0(cell_wrap_0 *pStruct);

extern void emxFreeStruct_cell_wrap_1(cell_wrap_1 *pStruct);

extern void emxFreeStruct_cell_wrap_2(cell_wrap_2 *pStruct);

extern void emxFreeStruct_cell_wrap_3(cell_wrap_3 *pStruct);

extern void emxFreeStruct_cell_wrap_4(cell_wrap_4 *pStruct);

extern void emxFreeStruct_cell_wrap_5(cell_wrap_5 *pStruct);

extern void emxFreeStruct_cell_wrap_6(cell_wrap_6 *pStruct);

extern void emxFreeStruct_cell_wrap_7(cell_wrap_7 *pStruct);

extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);

extern void emxFree_cell_wrap_0(emxArray_cell_wrap_0 **pEmxArray);

extern void emxFree_cell_wrap_1(emxArray_cell_wrap_1 **pEmxArray);

extern void emxFree_cell_wrap_2(emxArray_cell_wrap_2 **pEmxArray);

extern void emxFree_cell_wrap_3(emxArray_cell_wrap_3 **pEmxArray);

extern void emxFree_cell_wrap_4(emxArray_cell_wrap_4 **pEmxArray);

extern void emxFree_cell_wrap_5(emxArray_cell_wrap_5 **pEmxArray);

extern void emxFree_cell_wrap_6(emxArray_cell_wrap_6 **pEmxArray);

extern void emxFree_cell_wrap_7(emxArray_cell_wrap_7 **pEmxArray);

extern void emxFree_char_T(emxArray_char_T **pEmxArray);

extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);

extern void emxFree_real_T(emxArray_real_T **pEmxArray);

extern void emxInitStruct_cell_wrap_0(cell_wrap_0 *pStruct);

extern void emxInitStruct_cell_wrap_1(cell_wrap_1 *pStruct);

extern void emxInitStruct_cell_wrap_2(cell_wrap_2 *pStruct);

extern void emxInitStruct_cell_wrap_3(cell_wrap_3 *pStruct);

extern void emxInitStruct_cell_wrap_4(cell_wrap_4 *pStruct);

extern void emxInitStruct_cell_wrap_5(cell_wrap_5 *pStruct);

extern void emxInitStruct_cell_wrap_6(cell_wrap_6 *pStruct);

extern void emxInitStruct_cell_wrap_7(cell_wrap_7 *pStruct);

extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray,
                              int numDimensions);

extern void emxInit_cell_wrap_0(emxArray_cell_wrap_0 **pEmxArray);

extern void emxInit_cell_wrap_1(emxArray_cell_wrap_1 **pEmxArray);

extern void emxInit_cell_wrap_2(emxArray_cell_wrap_2 **pEmxArray);

extern void emxInit_cell_wrap_3(emxArray_cell_wrap_3 **pEmxArray);

extern void emxInit_cell_wrap_4(emxArray_cell_wrap_4 **pEmxArray);

extern void emxInit_cell_wrap_5(emxArray_cell_wrap_5 **pEmxArray);

extern void emxInit_cell_wrap_6(emxArray_cell_wrap_6 **pEmxArray,
                                int numDimensions);

extern void emxInit_cell_wrap_7(emxArray_cell_wrap_7 **pEmxArray);

extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);

extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);

extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

extern void emxTrim_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_2(emxArray_cell_wrap_2 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_3(emxArray_cell_wrap_3 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_4(emxArray_cell_wrap_4 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_5(emxArray_cell_wrap_5 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_6(emxArray_cell_wrap_6 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_7(emxArray_cell_wrap_7 *emxArray, int fromIndex,
                                int toIndex);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (mapTF_emxutil.h) */
