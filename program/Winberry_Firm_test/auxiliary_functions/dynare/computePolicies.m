function [vCapitalAdjust,vCapitalConstrained,vCutoff] = computePolicies(vCoefficients,vCapitalCoefficients,wage,mGrid,mPoly,aProdPrimePoly);

% Computes policy functions over arbitrary grid, mGrid
% 
% Inputs
% (1) vCoefficients: coefficients on value function
%	(2) vCapitalCoefficients: coefficients on capital accumulation policy, conditional on adjusting
%	(3) wage: wage
%	(4) mGrid: grid to compute policy on (size nGrid x 2)
%	(5) mPoly: polynomials over mGrid
%	(6) aProdPrimePoly: polynomials over future productivity
%
% Outputs
% (1) vCapitalAdjust: adjust capital decision
%	(2) vCapitalConstrained: constrained capital decision
%	(3) vCutoff: fixed cost threshold
% 
% Thomas Winberry, Feburary 14th, 2018

% Declare global variables used in this function
global ttheta nnu ddelta vStatePolySquared tolerance maxIterations nState ...
	aaUpper aaLower capitalMin capitalMax nShocks nProd nCapital bbeta ppsi2Capital dampening ...
	nProdFine nCapitalFine vShocksWeights ppsiCapital
	
% Get size of grid
[nGrid,~] 		= size(mGrid);


%----------------------------------------------------------------
% Capital accumulation policy
%----------------------------------------------------------------
	
% Conditional on adjusting
vCapitalAdjust = mPoly * vCapitalCoefficients;

% Conditional on not adjusting
vCapitalConstrained = vCapitalAdjust;
vCapitalConstrained(vCapitalAdjust > (1 - ddelta + aaUpper) * mGrid(:,2)) = ...
	(1 - ddelta + aaUpper) * mGrid(vCapitalAdjust > (1 - ddelta + aaUpper) * mGrid(:,2),2);
vCapitalConstrained(vCapitalAdjust < (1 - ddelta + aaLower) * mGrid(:,2)) = ...
	(1 - ddelta + aaLower) * mGrid(vCapitalAdjust < (1 - ddelta + aaLower) * mGrid(:,2),2);
	
	
%----------------------------------------------------------------
% Compute adjustment threshold
%----------------------------------------------------------------

% Load polynomials over next period's productivity
mProdPrimePoly 		= reshape(aProdPrimePoly,nShocks * nGrid,nProd);

% Polynomials over next period's capital, if adjust
vCapital 				= reshape(repmat(scaleDown(vCapitalAdjust,capitalMin,capitalMax)',[nShocks 1]),nShocks * nGrid,1);
mCapitalPrimePoly 	= computeChebyshev(nCapital,vCapital);
mCapitalAdjustPrimePoly = zeros(nShocks * nGrid,nState);
for iState = 1:nShocks * nGrid
	[vProd,vCapital] 							= ndgrid(mProdPrimePoly(iState,:),mCapitalPrimePoly(iState,:));
	mCapitalAdjustPrimePoly(iState,:) = reshape(vProd .* vCapital,nState,1);
end

% Polynomials over next period's capital, if don't adjust
vCapital 				= reshape(repmat(scaleDown(vCapitalConstrained,capitalMin,capitalMax)',[nShocks 1]),nShocks * nGrid,1);
mCapitalPrimePoly 	= computeChebyshev(nCapital,vCapital);
mCapitalConstrainedPrimePoly = zeros(nShocks * nGrid,nState);
for iState = 1:nShocks * nGrid
	[vProd,vCapital] 	= ndgrid(mProdPrimePoly(iState,:),mCapitalPrimePoly(iState,:));
	mCapitalConstrainedPrimePoly(iState,:) = reshape(vProd .* vCapital,nState,1);
end

% Compute expected value function, if adjust
mValuePrimeAdjust 		= reshape(mCapitalAdjustPrimePoly * vCoefficients,nShocks,nGrid);
vValueAdjust 					= (vShocksWeights' * mValuePrimeAdjust)';

% Compute expected value function, if don't adjust
mValuePrimeConstrained = reshape(mCapitalConstrainedPrimePoly * vCoefficients,nShocks,nGrid);
vValueConstrained 			= (vShocksWeights' * mValuePrimeConstrained)';

% Compute fixed cost threshold for adjustment
vCutoff = (1 / wage) * (-(vCapitalAdjust - vCapitalConstrained) + bbeta * (vValueAdjust - vValueConstrained));
	
% Ensure vCutoff stays on bounds of grid space
vCutoff = min(max(zeros(nGrid,1),vCutoff),ppsiCapital * ones(nGrid,1));
