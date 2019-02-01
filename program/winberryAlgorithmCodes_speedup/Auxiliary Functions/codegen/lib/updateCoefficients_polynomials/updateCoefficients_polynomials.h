/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * updateCoefficients_polynomials.h
 *
 * Code generation for function 'updateCoefficients_polynomials'
 *
 */

#ifndef UPDATECOEFFICIENTS_POLYNOMIALS_H
#define UPDATECOEFFICIENTS_POLYNOMIALS_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "updateCoefficients_polynomials_types.h"

/* Function Declarations */
extern void updateCoefficients_polynomials(const double mCoefficients[50],
  double bbeta, double ssigma, double aaBar, double mmu, double ttau, const
  double mEpsilonTransition[4], double nEpsilon, double nAssets, double nState,
  double assetsMin, double assetsMax, const double mEpsilonGrid[50], const
  double mAssetsGrid[50], const double vAssetsPoly[625], const double
  vAssetsPolySquared[25], const double mEpsilonPrimeGrid[100], double r, double
  w, emxArray_real_T *mCoefficientsNew);

#endif

/* End of code generation (updateCoefficients_polynomials.h) */
