/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_updateCoefficients_polynomials_api.c
 *
 * Code generation for function '_coder_updateCoefficients_polynomials_api'
 *
 */

/* Include files */
#include "tmwtypes.h"
#include "_coder_updateCoefficients_polynomials_api.h"
#include "_coder_updateCoefficients_polynomials_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131451U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "updateCoefficients_polynomials",    /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[50];
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *bbeta,
  const char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *mEpsilonTransition, const char_T *identifier))[4];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *mCoefficients, const char_T *identifier))[50];
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4];
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *vAssetsPoly, const char_T *identifier))[625];
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[625];
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *vAssetsPolySquared, const char_T *identifier))[25];
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[25];
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *mEpsilonPrimeGrid, const char_T *identifier))[100];
static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100];
static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[50];
static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4];
static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[625];
static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[25];
static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100];

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[50]
{
  real_T (*y)[50];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *bbeta,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(bbeta), &thisId);
  emlrtDestroyArray(&bbeta);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *mEpsilonTransition, const char_T *identifier))[4]
{
  real_T (*y)[4];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(mEpsilonTransition), &thisId);
  emlrtDestroyArray(&mEpsilonTransition);
  return y;
}
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *mCoefficients, const char_T *identifier))[50]
{
  real_T (*y)[50];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(mCoefficients), &thisId);
  emlrtDestroyArray(&mCoefficients);
  return y;
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 0, 0 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m0, (void *)&u->data[0]);
  emlrtSetDimensions((mxArray *)m0, u->size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex((void *)(*pEmxArray)->data);
    }

    emlrtFreeMex((void *)(*pEmxArray)->size);
    emlrtFreeMex((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)emlrtMallocMex(sizeof(emxArray_real_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2012b(sp, (void *)pEmxArray, (void (*)(void *))
      emxFree_real_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex((uint32_T)(sizeof(int32_T)
    * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[4]
{
  real_T (*y)[4];
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *vAssetsPoly, const char_T *identifier))[625]
{
  real_T (*y)[625];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(vAssetsPoly), &thisId);
  emlrtDestroyArray(&vAssetsPoly);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[625]
{
  real_T (*y)[625];
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *vAssetsPolySquared, const char_T *identifier))[25]
{
  real_T (*y)[25];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(vAssetsPolySquared), &thisId);
  emlrtDestroyArray(&vAssetsPolySquared);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[25]
{
  real_T (*y)[25];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *mEpsilonPrimeGrid, const char_T *identifier))[100]
{
  real_T (*y)[100];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(mEpsilonPrimeGrid), &thisId);
  emlrtDestroyArray(&mEpsilonPrimeGrid);
  return y;
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100]
{
  real_T (*y)[100];
  y = r_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[50]
{
  real_T (*ret)[50];
  static const int32_T dims[2] = { 2, 25 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[50])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[4]
{
  real_T (*ret)[4];
  static const int32_T dims[2] = { 2, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[625]
{
  real_T (*ret)[625];
  static const int32_T dims[2] = { 25, 25 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[625])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[25]
{
  real_T (*ret)[25];
  static const int32_T dims[1] = { 25 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[25])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100]
{
  real_T (*ret)[100];
  static const int32_T dims[2] = { 2, 50 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[100])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void updateCoefficients_polynomials_api(const mxArray * const prhs[19], const
  mxArray *plhs[1])
{
  emxArray_real_T *mCoefficientsNew;
  real_T (*mCoefficients)[50];
  real_T bbeta;
  real_T ssigma;
  real_T aaBar;
  real_T mmu;
  real_T ttau;
  real_T (*mEpsilonTransition)[4];
  real_T nEpsilon;
  real_T nAssets;
  real_T nState;
  real_T assetsMin;
  real_T assetsMax;
  real_T (*mEpsilonGrid)[50];
  real_T (*mAssetsGrid)[50];
  real_T (*vAssetsPoly)[625];
  real_T (*vAssetsPolySquared)[25];
  real_T (*mEpsilonPrimeGrid)[100];
  real_T r;
  real_T w;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &mCoefficientsNew, 2, true);

  /* Marshall function inputs */
  mCoefficients = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "mCoefficients");
  bbeta = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "bbeta");
  ssigma = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "ssigma");
  aaBar = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "aaBar");
  mmu = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "mmu");
  ttau = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "ttau");
  mEpsilonTransition = e_emlrt_marshallIn(&st, emlrtAlias(prhs[6]),
    "mEpsilonTransition");
  nEpsilon = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "nEpsilon");
  nAssets = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "nAssets");
  nState = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "nState");
  assetsMin = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "assetsMin");
  assetsMax = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "assetsMax");
  mEpsilonGrid = emlrt_marshallIn(&st, emlrtAlias(prhs[12]), "mEpsilonGrid");
  mAssetsGrid = emlrt_marshallIn(&st, emlrtAlias(prhs[13]), "mAssetsGrid");
  vAssetsPoly = g_emlrt_marshallIn(&st, emlrtAlias(prhs[14]), "vAssetsPoly");
  vAssetsPolySquared = i_emlrt_marshallIn(&st, emlrtAlias(prhs[15]),
    "vAssetsPolySquared");
  mEpsilonPrimeGrid = k_emlrt_marshallIn(&st, emlrtAlias(prhs[16]),
    "mEpsilonPrimeGrid");
  r = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[17]), "r");
  w = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[18]), "w");

  /* Invoke the target function */
  updateCoefficients_polynomials(*mCoefficients, bbeta, ssigma, aaBar, mmu, ttau,
    *mEpsilonTransition, nEpsilon, nAssets, nState, assetsMin, assetsMax,
    *mEpsilonGrid, *mAssetsGrid, *vAssetsPoly, *vAssetsPolySquared,
    *mEpsilonPrimeGrid, r, w, mCoefficientsNew);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(mCoefficientsNew);
  mCoefficientsNew->canFreeData = false;
  emxFree_real_T(&mCoefficientsNew);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

void updateCoefficients_polynomials_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  updateCoefficients_polynomials_xil_terminate();
}

void updateCoefficients_polynomials_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void updateCoefficients_polynomials_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_updateCoefficients_polynomials_api.c) */
