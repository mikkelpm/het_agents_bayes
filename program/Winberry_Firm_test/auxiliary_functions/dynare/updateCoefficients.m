function [vCoefficientsNew,vCapitalAdjust] = updateCoefficients(vCoefficients)

% Updates coefficients approximating value function in steady state
% 
% Inputs
%	(1) vCoefficients: coefficients of value function
%
% Outputs 
%	(1) vCoefficientsNew: updated coefficients of value function
%	(2) vCapitalAdjust: capital decision conditional on adjusting, given vCoefficients
%
% Thomas Winberry, Feburary 14th, 2018

% Declare global variables used in this function
global ddelta aaUpper aaLower ppsiCapital ...
	bbeta nProd nCapital nState capitalMin capitalMax ...
	nShocks acc ...
	vShocksWeights mStateGrid mStatePoly vStatePolySquared aProdPrimePoly ...
	wage vProfitGrid 
	
%---------------------------------------------------------------
% Compute investment decision conditional on adjustment or not
%---------------------------------------------------------------

% Compute decision if adjust (only do it once per productivity value since its independent of capital)
vCapitalAdjust 		= zeros(nState,1);
options 					= optimset('TolX',1e-12);
vCapitalAdjustProd 	= zeros(nProd,1);

for iEpsilon = 1:nProd
	f 												= @(kprime) capitalResid(kprime,squeeze(aProdPrimePoly(:,iEpsilon,:)),vCoefficients);
	vCapitalAdjustEpsilon(iEpsilon,1) 	= fminbnd(f,capitalMin,capitalMax,options);	
end

% Expand along entire grid, including capital
vCapitalAdjust 								= reshape(vCapitalAdjustEpsilon * ones(1,nCapital),nState,1);

% Policy function if do not adjust
vCapitalConstrained 						= vCapitalAdjust;
vCapitalConstrained(vCapitalAdjust > (1 - ddelta + aaUpper) * mStateGrid(:,2)) = ...
	(1 - ddelta + aaUpper) * mStateGrid(vCapitalAdjust > (1 - ddelta + aaUpper) * mStateGrid(:,2),2);
vCapitalConstrained(vCapitalAdjust < (1 - ddelta + aaLower) * mStateGrid(:,2)) = ...
	(1 - ddelta + aaLower) * mStateGrid(vCapitalAdjust < (1 - ddelta + aaLower) * mStateGrid(:,2),2);
	
	
%---------------------------------------------------------------
% Compute polynomials of next period's value function, 
% conditional on investment decision
%---------------------------------------------------------------

% Load polynomials over next period's productivity
mProdPrimePoly 		= reshape(aProdPrimePoly,nShocks * nState,nProd);

% Polynomials over next period's capital, if adjust
vCapital 				= reshape(repmat(scaleDown(vCapitalAdjust,capitalMin,capitalMax)',[nShocks 1]),nShocks * nState,1);
mCapitalPrimePoly	= computeChebyshev(nCapital,vCapital);
mCapitalAdjustPrimePoly = zeros(nShocks * nState,nState);
for iState = 1:nShocks * nState
	[vProd,vCapital] 							= ndgrid(mProdPrimePoly(iState,:),mCapitalPrimePoly(iState,:));
	mCapitalAdjustPrimePoly(iState,:) = reshape(vProd .* vCapital,nState,1);
end

% Polynomials over next period's capital, if don't adjust
vCapital 								= reshape(repmat(scaleDown(vCapitalConstrained,capitalMin,capitalMax)',[nShocks 1]),nShocks * nState,1);
mCapitalPrimePoly 					= computeChebyshev(nCapital,vCapital);
mCapitalConstrainedPrimePoly 	= zeros(nShocks * nState,nState);
for iState = 1:nShocks * nState
	[vProd,vCapital] 					= ndgrid(mProdPrimePoly(iState,:),mCapitalPrimePoly(iState,:));
	mCapitalConstrainedPrimePoly(iState,:) = reshape(vProd .* vCapital,nState,1);
end


%---------------------------------------------------------------
% Update value function "acc" times using these policy functions
%---------------------------------------------------------------

iterations = 1;
while iterations <= acc

	% Compute expected value function, if adjust
	mValuePrimeAdjust 		= reshape(mCapitalAdjustPrimePoly * vCoefficients,nShocks,nState);
	vValueAdjust 					= (vShocksWeights' * mValuePrimeAdjust)';
	
	% Compute expected value function, if don't adjust
	mValuePrimeConstrained 	= reshape(mCapitalConstrainedPrimePoly * vCoefficients,nShocks,nState);
	vValueConstrained			= (vShocksWeights' * mValuePrimeConstrained)';
	
	% Compute fixed cost threshold for adjustment
	vCutoff 							= (1 / wage) * (-(vCapitalAdjust - vCapitalConstrained) + bbeta * (vValueAdjust - vValueConstrained));
	vCutoff 							= min(max(zeros(nState,1),vCutoff),ppsiCapital * ones(nState,1));		% enforce grid bounds
	
	% Compute right hand side of Bellman operator
	vNewGrid 						= vProfitGrid + (vCutoff ./ ppsiCapital) .* (-(vCapitalAdjust - (1 - ddelta) * mStateGrid(:,2)) - ...
											wage * (vCutoff / 2) + bbeta * vValueAdjust) + (1 - (vCutoff ./ ppsiCapital)) .* (-(vCapitalConstrained - ...
											(1 - ddelta) * mStateGrid(:,2)) + bbeta * vValueConstrained);
	
	% Compute new coefficients
	vCoefficients = sum(mStatePoly' .* (ones(nState,1) * vNewGrid'),2);
	vCoefficients = vCoefficients ./ vStatePolySquared;
	
	% Update iteration counter
	iterations = iterations + 1;

end


%---------------------------------------------------------------
% Output of the function
%---------------------------------------------------------------

vCoefficientsNew = vCoefficients;
