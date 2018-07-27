function value = capitalResid(kprime,mProdPrimePoly,vCoefficients)

% Objective function for optimization problem of capital accumulation, conditional on adjusting
% 
% Inputs
% (1) kprime: candidate choice of future capital
%	(2) mProdPrimePoly: polynomial over future shocks, conditional on realization of current shock
%	(3) vCoefficients: coefficients of value function
%
% Outputs
%  (1) value: negative of value associated with choice kprime (since using minimizer rather than maximizer)
% 
% Thomas Winberry, February 14th, 2018

% Declare global variables
global nCapital capitalMin capitalMax nShocks nState bbeta vShocksWeights ddelta

% Compute polynomial over capital choice
mCapitalPrimePoly 	= computeChebyshev(nCapital,scaleDown(kprime,capitalMin,capitalMax)*ones(nShocks,1));

% Compute tensor product of polynomials
mPolyPrime			= zeros(nShocks,nState);
for iShock = 1:nShocks
	[vProd,vCapital] 			= ndgrid(mProdPrimePoly(iShock,:),mCapitalPrimePoly(iShock,:));
	mPolyPrime(iShock,:) 	= reshape(vProd .* vCapital,nState,1);
end

% Compute next period's value function
vValuePrime 			= mPolyPrime * vCoefficients;

% Compute expectation
expectedValuePrime = vShocksWeights' * vValuePrime;

% Compute objective function
value 					= -(-kprime + bbeta * expectedValuePrime);