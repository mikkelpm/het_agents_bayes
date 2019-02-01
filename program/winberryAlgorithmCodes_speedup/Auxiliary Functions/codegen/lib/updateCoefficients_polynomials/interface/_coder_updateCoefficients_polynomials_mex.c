/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_updateCoefficients_polynomials_mex.c
 *
 * Code generation for function '_coder_updateCoefficients_polynomials_mex'
 *
 */

/* Include files */
#include "_coder_updateCoefficients_polynomials_api.h"
#include "_coder_updateCoefficients_polynomials_mex.h"

/* Function Declarations */
static void c_updateCoefficients_polynomial(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[19]);

/* Function Definitions */
static void c_updateCoefficients_polynomial(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[19])
{
  const mxArray *inputs[19];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 19) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 19, 4,
                        30, "updateCoefficients_polynomials");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 30,
                        "updateCoefficients_polynomials");
  }

  /* Temporary copy for mex inputs. */
  if (0 <= nrhs - 1) {
    memcpy((void *)&inputs[0], (void *)&prhs[0], (uint32_T)(nrhs * (int32_T)
            sizeof(const mxArray *)));
  }

  /* Call the function. */
  updateCoefficients_polynomials_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  updateCoefficients_polynomials_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(updateCoefficients_polynomials_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  updateCoefficients_polynomials_initialize();

  /* Dispatch the entry-point. */
  c_updateCoefficients_polynomial(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_updateCoefficients_polynomials_mex.c) */
