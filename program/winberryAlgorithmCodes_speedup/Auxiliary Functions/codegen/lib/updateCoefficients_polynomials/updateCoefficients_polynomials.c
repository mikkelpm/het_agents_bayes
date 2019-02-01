/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * updateCoefficients_polynomials.c
 *
 * Code generation for function 'updateCoefficients_polynomials'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "updateCoefficients_polynomials.h"
#include "updateCoefficients_polynomials_emxutil.h"
#include "repmat.h"
#include "power.h"
#include "exp.h"
#include "updateCoefficients_polynomials_rtwutil.h"

/* Function Definitions */
void updateCoefficients_polynomials(const double mCoefficients[50], double bbeta,
  double ssigma, double aaBar, double mmu, double ttau, const double
  mEpsilonTransition[4], double nEpsilon, double nAssets, double nState, double
  assetsMin, double assetsMax, const double mEpsilonGrid[50], const double
  mAssetsGrid[50], const double vAssetsPoly[625], const double
  vAssetsPolySquared[25], const double mEpsilonPrimeGrid[100], double r, double
  w, emxArray_real_T *mCoefficientsNew)
{
  int ia;
  double mAssetsPrime[50];
  int ar;
  double mAssetsPrimeStar[50];
  int k;
  int nx;
  emxArray_real_T *mAssetsPrimeGrid;
  double b_mAssetsPrimeStar;
  double b_nEpsilon[2];
  int iv0[2];
  double b_mAssetsPrime;
  emxArray_real_T *mPolyAssetsPrime;
  int loop_ub;
  emxArray_int32_T *r0;
  emxArray_real_T *b_mPolyAssetsPrime;
  emxArray_real_T *c_mPolyAssetsPrime;
  emxArray_real_T *mConditionalExpectationPrime;
  emxArray_real_T *b;
  unsigned int unnamed_idx_1;
  emxArray_real_T *z1;
  int br;
  int ic;
  double mAssetsPrimePrimeStar[100];
  double mAssetsPrimePrimeGrid[100];
  double b_z1[100];
  emxArray_real_T *aConditionalExpectationTilde;
  emxArray_real_T *x;
  emxArray_real_T *b_b;
  double a[625];
  double vCoefficients[25];
  double b_vCoefficients[25];

  /*  Updates polynomial coefficients approximating the conditional expectation function in steady state */
  /*   */
  /*  Inputs */
  /*    (1) mCoefficients: nEpsilon x nAssets matrix, storing previous iteration's coefficients */
  /*  */
  /*  Outputs */
  /*    (1) mCoefficientsNew: nEpsilon x nAssets matrix, storing updated coefficients */
  /*   */
  /*  Thomas Winberry, January 19, 2016 */
  /*  Declare global variables */
  /* --------------------------------------------------------------- */
  /*  Compute current period's savings policy function */
  /* --------------------------------------------------------------- */
  /*  Compute conditional expectation */
  /*  Compute target saving */
  for (ia = 0; ia < 2; ia++) {
    for (ar = 0; ar < 25; ar++) {
      mAssetsPrime[ia + (ar << 1)] = 0.0;
      for (nx = 0; nx < 25; nx++) {
        mAssetsPrime[ia + (ar << 1)] += mCoefficients[ia + (nx << 1)] *
          vAssetsPoly[ar + 25 * nx];
      }
    }
  }

  b_exp(mAssetsPrime);
  power(mAssetsPrime, -1.0 / ssigma, mAssetsPrimeStar);

  /*  Compute actual saving */
  for (k = 0; k < 50; k++) {
    b_mAssetsPrimeStar = (w * (mmu * (1.0 - mEpsilonGrid[k]) + (1.0 - ttau) *
      mEpsilonGrid[k]) + (1.0 + r) * mAssetsGrid[k]) - mAssetsPrimeStar[k];
    if ((b_mAssetsPrimeStar > aaBar) || rtIsNaN(aaBar)) {
      mAssetsPrime[k] = b_mAssetsPrimeStar;
    } else {
      mAssetsPrime[k] = aaBar;
    }
  }

  emxInit_real_T(&mAssetsPrimeGrid, 2);
  b_nEpsilon[0] = nEpsilon;
  b_nEpsilon[1] = 1.0;
  iv0[0] = 1;
  iv0[1] = (int)nState;
  repmat(mAssetsPrime, iv0, b_nEpsilon, mAssetsPrimeGrid);

  /*  Compute next period's polynomials */
  b_mAssetsPrimeStar = assetsMax - assetsMin;
  for (k = 0; k < 50; k++) {
    b_mAssetsPrime = 2.0 * ((mAssetsPrime[k] - assetsMin) / b_mAssetsPrimeStar)
      - 1.0;
    if (!(b_mAssetsPrime > -1.0)) {
      b_mAssetsPrime = -1.0;
    }

    mAssetsPrimeStar[k] = b_mAssetsPrime;
  }

  for (k = 0; k < 50; k++) {
    b_mAssetsPrimeStar = mAssetsPrimeStar[k];
    if (!(b_mAssetsPrimeStar < 1.0)) {
      b_mAssetsPrimeStar = 1.0;
    }

    mAssetsPrime[k] = b_mAssetsPrimeStar;
  }

  emxInit_real_T(&mPolyAssetsPrime, 2);
  ia = mPolyAssetsPrime->size[0] * mPolyAssetsPrime->size[1];
  mPolyAssetsPrime->size[0] = (int)nState;
  mPolyAssetsPrime->size[1] = (int)nAssets;
  emxEnsureCapacity_real_T(mPolyAssetsPrime, ia);
  loop_ub = (int)nState * (int)nAssets;
  for (ia = 0; ia < loop_ub; ia++) {
    mPolyAssetsPrime->data[ia] = 1.0;
  }

  emxInit_int32_T(&r0, 1);
  ia = r0->size[0];
  r0->size[0] = (int)nState;
  emxEnsureCapacity_int32_T(r0, ia);
  loop_ub = (int)nState;
  for (ia = 0; ia < loop_ub; ia++) {
    r0->data[ia] = ia;
  }

  nx = r0->size[0];
  for (ia = 0; ia < nx; ia++) {
    mPolyAssetsPrime->data[r0->data[ia] + mPolyAssetsPrime->size[0]] =
      mAssetsPrime[ia];
  }

  k = 0;
  emxInit_real_T2(&b_mPolyAssetsPrime, 1);
  emxInit_real_T2(&c_mPolyAssetsPrime, 1);
  while (k <= (int)(nAssets + -2.0) - 1) {
    loop_ub = mPolyAssetsPrime->size[0];
    ia = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(r0, ia);
    for (ia = 0; ia < loop_ub; ia++) {
      r0->data[ia] = ia;
    }

    nx = mPolyAssetsPrime->size[0];
    ia = b_mPolyAssetsPrime->size[0];
    b_mPolyAssetsPrime->size[0] = nx;
    emxEnsureCapacity_real_T1(b_mPolyAssetsPrime, ia);
    for (ia = 0; ia < nx; ia++) {
      b_mPolyAssetsPrime->data[ia] = mPolyAssetsPrime->data[ia +
        mPolyAssetsPrime->size[0] * ((int)((3.0 + (double)k) - 1.0) - 1)];
    }

    nx = mPolyAssetsPrime->size[0];
    ia = c_mPolyAssetsPrime->size[0];
    c_mPolyAssetsPrime->size[0] = nx;
    emxEnsureCapacity_real_T1(c_mPolyAssetsPrime, ia);
    for (ia = 0; ia < nx; ia++) {
      c_mPolyAssetsPrime->data[ia] = mPolyAssetsPrime->data[ia +
        mPolyAssetsPrime->size[0] * ((int)((3.0 + (double)k) - 2.0) - 1)];
    }

    for (ia = 0; ia < 50; ia++) {
      mAssetsPrimeStar[ia] = 2.0 * mAssetsPrime[ia] * b_mPolyAssetsPrime->
        data[ia] - c_mPolyAssetsPrime->data[ia];
    }

    nx = r0->size[0];
    for (ia = 0; ia < nx; ia++) {
      mPolyAssetsPrime->data[r0->data[ia] + mPolyAssetsPrime->size[0] * (k + 2)]
        = mAssetsPrimeStar[ia];
    }

    k++;
  }

  emxFree_real_T(&c_mPolyAssetsPrime);
  emxFree_real_T(&b_mPolyAssetsPrime);
  emxInit_real_T(&mConditionalExpectationPrime, 2);
  emxInit_real_T(&b, 2);

  /* --------------------------------------------------------------- */
  /*  Compute next period's savings policy function */
  /* --------------------------------------------------------------- */
  /*  Compute conditional expectation */
  ia = b->size[0] * b->size[1];
  b->size[0] = mPolyAssetsPrime->size[1];
  b->size[1] = mPolyAssetsPrime->size[0];
  emxEnsureCapacity_real_T(b, ia);
  loop_ub = mPolyAssetsPrime->size[0];
  for (ia = 0; ia < loop_ub; ia++) {
    k = mPolyAssetsPrime->size[1];
    for (ar = 0; ar < k; ar++) {
      b->data[ar + b->size[0] * ia] = mPolyAssetsPrime->data[ia +
        mPolyAssetsPrime->size[0] * ar];
    }
  }

  if (b->size[0] == 1) {
    ia = mConditionalExpectationPrime->size[0] *
      mConditionalExpectationPrime->size[1];
    mConditionalExpectationPrime->size[0] = 2;
    mConditionalExpectationPrime->size[1] = b->size[1];
    emxEnsureCapacity_real_T(mConditionalExpectationPrime, ia);
    for (ia = 0; ia < 2; ia++) {
      loop_ub = b->size[1];
      for (ar = 0; ar < loop_ub; ar++) {
        mConditionalExpectationPrime->data[ia +
          mConditionalExpectationPrime->size[0] * ar] = 0.0;
        for (nx = 0; nx < 25; nx++) {
          mConditionalExpectationPrime->data[ia +
            mConditionalExpectationPrime->size[0] * ar] += mCoefficients[ia +
            (nx << 1)] * b->data[nx + b->size[0] * ar];
        }
      }
    }
  } else {
    unnamed_idx_1 = (unsigned int)b->size[1];
    ia = mConditionalExpectationPrime->size[0] *
      mConditionalExpectationPrime->size[1];
    mConditionalExpectationPrime->size[1] = (int)unnamed_idx_1;
    mConditionalExpectationPrime->size[0] = 2;
    emxEnsureCapacity_real_T(mConditionalExpectationPrime, ia);
    loop_ub = mConditionalExpectationPrime->size[1];
    for (ia = 0; ia < loop_ub; ia++) {
      for (ar = 0; ar < 2; ar++) {
        mConditionalExpectationPrime->data[ar +
          mConditionalExpectationPrime->size[0] * ia] = 0.0;
      }
    }

    if (b->size[1] != 0) {
      nx = (b->size[1] - 1) << 1;
      for (k = 0; k <= nx; k += 2) {
        for (ic = k; ic + 1 <= k + 2; ic++) {
          mConditionalExpectationPrime->data[ic] = 0.0;
        }
      }

      br = 0;
      for (k = 0; k <= nx; k += 2) {
        ar = 0;
        for (loop_ub = br; loop_ub + 1 <= br + 25; loop_ub++) {
          if (b->data[loop_ub] != 0.0) {
            ia = ar;
            for (ic = k; ic + 1 <= k + 2; ic++) {
              ia++;
              mConditionalExpectationPrime->data[ic] += b->data[loop_ub] *
                mCoefficients[ia - 1];
            }
          }

          ar += 2;
        }

        br += 25;
      }
    }
  }

  emxFree_real_T(&b);
  nx = mConditionalExpectationPrime->size[1] << 1;
  for (k = 0; k + 1 <= nx; k++) {
    mConditionalExpectationPrime->data[k] = exp
      (mConditionalExpectationPrime->data[k]);
  }

  emxInit_real_T(&z1, 2);

  /*  Compute target saving */
  b_mAssetsPrimeStar = -1.0 / ssigma;
  ia = z1->size[0] * z1->size[1];
  z1->size[0] = 2;
  z1->size[1] = mConditionalExpectationPrime->size[1];
  emxEnsureCapacity_real_T(z1, ia);
  unnamed_idx_1 = (unsigned int)mConditionalExpectationPrime->size[1];
  nx = (int)unnamed_idx_1 << 1;
  for (k = 0; k + 1 <= nx; k++) {
    z1->data[k] = rt_powd_snf(mConditionalExpectationPrime->data[k],
      b_mAssetsPrimeStar);
  }

  emxFree_real_T(&mConditionalExpectationPrime);
  for (ia = 0; ia < 100; ia++) {
    mAssetsPrimePrimeStar[ia] = (w * (mmu * (1.0 - mEpsilonPrimeGrid[ia]) + (1.0
      - ttau) * mEpsilonPrimeGrid[ia]) + (1.0 + r) * mAssetsPrimeGrid->data[ia])
      - z1->data[ia];
  }

  emxFree_real_T(&z1);

  /*  Compute actual savings */
  /*  mAssetsPrimePrimeGrid = max(mAssetsPrimePrimeStar,aaBar*ones(nEpsilon,nEpsilon*nAssets)); */
  for (k = 0; k < 100; k++) {
    if ((mAssetsPrimePrimeStar[k] > aaBar) || rtIsNaN(aaBar)) {
      mAssetsPrimePrimeGrid[k] = mAssetsPrimePrimeStar[k];
    } else {
      mAssetsPrimePrimeGrid[k] = aaBar;
    }
  }

  /* --------------------------------------------------------------- */
  /*  Update conditional expectation function */
  /* --------------------------------------------------------------- */
  /*  Compute new conditional expectation function */
  for (ia = 0; ia < 100; ia++) {
    b_mAssetsPrimeStar = (w * (mmu * (1.0 - mEpsilonPrimeGrid[ia]) + (1.0 - ttau)
      * mEpsilonPrimeGrid[ia]) + (1.0 + r) * mAssetsPrimeGrid->data[ia]) -
      mAssetsPrimePrimeGrid[ia];
    mAssetsPrimePrimeGrid[ia] = b_mAssetsPrimeStar;
  }

  emxFree_real_T(&mAssetsPrimeGrid);
  for (k = 0; k < 100; k++) {
    b_z1[k] = rt_powd_snf(mAssetsPrimePrimeGrid[k], -ssigma);
  }

  for (ia = 0; ia < 2; ia++) {
    for (ar = 0; ar < 50; ar++) {
      mAssetsPrimePrimeStar[ia + (ar << 1)] = 0.0;
      for (nx = 0; nx < 2; nx++) {
        mAssetsPrimePrimeStar[ia + (ar << 1)] += bbeta * mEpsilonTransition[ia +
          (nx << 1)] * ((1.0 + r) * b_z1[nx + (ar << 1)]);
      }
    }
  }

  emxInit_real_T1(&aConditionalExpectationTilde, 3);
  ia = aConditionalExpectationTilde->size[0] *
    aConditionalExpectationTilde->size[1] * aConditionalExpectationTilde->size[2];
  aConditionalExpectationTilde->size[0] = (int)nEpsilon;
  aConditionalExpectationTilde->size[1] = (int)nEpsilon;
  aConditionalExpectationTilde->size[2] = (int)nAssets;
  emxEnsureCapacity_real_T2(aConditionalExpectationTilde, ia);
  loop_ub = (int)nAssets;
  for (ia = 0; ia < loop_ub; ia++) {
    k = (int)nEpsilon;
    for (ar = 0; ar < k; ar++) {
      br = (int)nEpsilon;
      for (nx = 0; nx < br; nx++) {
        aConditionalExpectationTilde->data[(nx +
          aConditionalExpectationTilde->size[0] * ar) +
          aConditionalExpectationTilde->size[0] *
          aConditionalExpectationTilde->size[1] * ia] = mAssetsPrimePrimeStar
          [(nx + (int)nEpsilon * ar) + (int)nEpsilon * (int)nEpsilon * ia];
      }
    }
  }

  /*  Extract the relevant entries */
  ia = mPolyAssetsPrime->size[0] * mPolyAssetsPrime->size[1];
  mPolyAssetsPrime->size[0] = (int)nEpsilon;
  mPolyAssetsPrime->size[1] = (int)nAssets;
  emxEnsureCapacity_real_T(mPolyAssetsPrime, ia);
  for (br = 0; br < (int)nEpsilon; br++) {
    loop_ub = mPolyAssetsPrime->size[1];
    ia = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(r0, ia);
    for (ia = 0; ia < loop_ub; ia++) {
      r0->data[ia] = ia;
    }

    loop_ub = aConditionalExpectationTilde->size[2];
    for (ia = 0; ia < loop_ub; ia++) {
      mAssetsPrimePrimeStar[ia] = aConditionalExpectationTilde->data[(br +
        aConditionalExpectationTilde->size[0] * br) +
        aConditionalExpectationTilde->size[0] *
        aConditionalExpectationTilde->size[1] * ia];
    }

    nx = r0->size[0];
    for (ia = 0; ia < nx; ia++) {
      mPolyAssetsPrime->data[br + mPolyAssetsPrime->size[0] * r0->data[ia]] =
        mAssetsPrimePrimeStar[ia];
    }
  }

  emxFree_real_T(&aConditionalExpectationTilde);

  /*  Update the coefficients */
  ia = mCoefficientsNew->size[0] * mCoefficientsNew->size[1];
  mCoefficientsNew->size[0] = (int)nEpsilon;
  mCoefficientsNew->size[1] = (int)nAssets;
  emxEnsureCapacity_real_T(mCoefficientsNew, ia);
  br = 0;
  emxInit_real_T(&x, 2);
  emxInit_real_T2(&b_b, 1);
  if (0 <= (int)nEpsilon - 1) {
    for (ia = 0; ia < 25; ia++) {
      for (ar = 0; ar < 25; ar++) {
        a[ar + 25 * ia] = vAssetsPoly[ia + 25 * ar];
      }
    }
  }

  while (br <= (int)nEpsilon - 1) {
    /*  	vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * log(mConditionalExpectation(iEpsilon,:))),2); */
    loop_ub = mPolyAssetsPrime->size[1];
    ia = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = loop_ub;
    emxEnsureCapacity_real_T(x, ia);
    for (ia = 0; ia < loop_ub; ia++) {
      x->data[x->size[0] * ia] = mPolyAssetsPrime->data[br +
        mPolyAssetsPrime->size[0] * ia];
    }

    ia = mPolyAssetsPrime->size[1];
    for (k = 0; k + 1 <= ia; k++) {
      x->data[k] = log(x->data[k]);
    }

    ia = b_b->size[0];
    b_b->size[0] = x->size[1];
    emxEnsureCapacity_real_T1(b_b, ia);
    loop_ub = x->size[1];
    for (ia = 0; ia < loop_ub; ia++) {
      b_b->data[ia] = x->data[x->size[0] * ia];
    }

    if (b_b->size[0] == 1) {
      for (ia = 0; ia < 25; ia++) {
        vCoefficients[ia] = 0.0;
        for (ar = 0; ar < 25; ar++) {
          b_mAssetsPrimeStar = vCoefficients[ia] + a[ia + 25 * ar] * b_b->
            data[ar];
          vCoefficients[ia] = b_mAssetsPrimeStar;
        }
      }
    } else if (25 == b_b->size[0]) {
      for (ia = 0; ia < 25; ia++) {
        vCoefficients[ia] = 0.0;
        for (ar = 0; ar < 25; ar++) {
          b_mAssetsPrimeStar = vCoefficients[ia] + a[ia + 25 * ar] * b_b->
            data[ar];
          vCoefficients[ia] = b_mAssetsPrimeStar;
        }
      }
    } else {
      memset(&vCoefficients[0], 0, 25U * sizeof(double));
      ar = 0;
      for (loop_ub = 0; loop_ub < 25; loop_ub++) {
        if (b_b->data[loop_ub] != 0.0) {
          ia = ar;
          for (ic = 0; ic < 25; ic++) {
            ia++;
            b_mAssetsPrimeStar = vCoefficients[ic] + b_b->data[loop_ub] * a[ia -
              1];
            vCoefficients[ic] = b_mAssetsPrimeStar;
          }
        }

        ar += 25;
      }
    }

    loop_ub = mCoefficientsNew->size[1];
    ia = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(r0, ia);
    for (ia = 0; ia < loop_ub; ia++) {
      r0->data[ia] = ia;
    }

    for (ia = 0; ia < 25; ia++) {
      b_vCoefficients[ia] = vCoefficients[ia] / vAssetsPolySquared[ia];
    }

    nx = r0->size[0];
    for (ia = 0; ia < nx; ia++) {
      mCoefficientsNew->data[br + mCoefficientsNew->size[0] * r0->data[ia]] =
        b_vCoefficients[ia];
    }

    br++;
  }

  emxFree_real_T(&b_b);
  emxFree_real_T(&x);
  emxFree_int32_T(&r0);
  emxFree_real_T(&mPolyAssetsPrime);
}

/* End of code generation (updateCoefficients_polynomials.c) */
