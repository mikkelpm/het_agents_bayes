/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_updateCoefficients_polynomials_api.h
 *
 * Code generation for function '_coder_updateCoefficients_polynomials_api'
 *
 */

#ifndef _CODER_UPDATECOEFFICIENTS_POLYNOMIALS_API_H
#define _CODER_UPDATECOEFFICIENTS_POLYNOMIALS_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_updateCoefficients_polynomials_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void updateCoefficients_polynomials(real_T mCoefficients[50], real_T
  bbeta, real_T ssigma, real_T aaBar, real_T mmu, real_T ttau, real_T
  mEpsilonTransition[4], real_T nEpsilon, real_T nAssets, real_T nState, real_T
  assetsMin, real_T assetsMax, real_T mEpsilonGrid[50], real_T mAssetsGrid[50],
  real_T vAssetsPoly[625], real_T vAssetsPolySquared[25], real_T
  mEpsilonPrimeGrid[100], real_T r, real_T w, emxArray_real_T *mCoefficientsNew);
extern void updateCoefficients_polynomials_api(const mxArray * const prhs[19],
  const mxArray *plhs[1]);
extern void updateCoefficients_polynomials_atexit(void);
extern void updateCoefficients_polynomials_initialize(void);
extern void updateCoefficients_polynomials_terminate(void);
extern void updateCoefficients_polynomials_xil_terminate(void);

#endif

/* End of code generation (_coder_updateCoefficients_polynomials_api.h) */
