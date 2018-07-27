% Compute polynomials for approximating various objects
%
% Thomas Winberry, February 14th, 2018

%---------------------------------------------------------------
% Polynomials for value function and capital accumulation conditional on adjusting
%		(Chebyshev polynomials evaluated over collocation nodes from computeGrids.m)
%---------------------------------------------------------------

global mStatePoly vStatePolySquared aProdPrimePoly

% Create one-dimensional polynomials
mProdPoly 		= computeChebyshev(nProd,mStateGridZeros(:,1));
mCapitalPoly 	= computeChebyshev(nCapital,mStateGridZeros(:,2));

% Create tensor product in polynomials
mStatePoly 	= zeros(nState,nState);	
for iState = 1:nState
	[vProd,vCapital] 			= ndgrid(mProdPoly(iState,:),mCapitalPoly(iState,:));
	mStatePoly(iState,:) 	= reshape(vProd .* vCapital,nState,1);
end
clear mProdPoly mCapitalPoly vProd vCapital

% Compute squared terms for interpolation formulas
mProdPoly 				= computeChebyshev(nProd,vProdZeros);
mCapitalPoly 			= computeChebyshev(nCapital,vCapitalZeros);
vProd 					= sum(mProdPoly .^ 2)';
vCapital 				= sum(mCapitalPoly .^ 2)';
vStatePolySquared 	= reshape(vProd * vCapital',nState,1);
clear mProdPoly mCapitalPoly vProd vCapital

% Compute polynomials over future productivity shocks (for computing expected value function)
aProdPrimePoly 		= reshape(computeChebyshev(nProd,reshape(mProdPrimeZeros,nShocks * nState,1)),nShocks,nState,nProd);


%---------------------------------------------------------------
% Chebyshev polynomials over fine grid
%		(useful for plotting and computing law of motion for histogram)
%---------------------------------------------------------------

global mFinePoly aProdPrimeFinePoly

% Create one-dimensional polynomials
mProdFinePoly 		= computeChebyshev(nProd,mFineGridZeros(:,1));
mCapitalFinePoly 	= computeChebyshev(nCapital,mFineGridZeros(:,2));

% Compute tensor product of polynomials
mFinePoly 				= zeros(nStateFine,nState);
for iState = 1:nStateFine
	[vProd,vCapital] 		= ndgrid(mProdFinePoly(iState,:),mCapitalFinePoly(iState,:));
	mFinePoly(iState,:) 	= reshape(vProd .* vCapital,nState,1);
end
clear vProd vCapital

% Compute polynomials over future shocks
aProdPrimeFinePoly 	= reshape(computeChebyshev(nProd,reshape(mProdPrimeFineZeros,nShocks * nStateFine,1)),nShocks,...
									nStateFine,nProd);
clear mProdFinePoly mCapitalFinePoly


%---------------------------------------------------------------
% Chebyshev polynomials over quadrature grid 
%		(useful for compting law of motion for parametric family)
%---------------------------------------------------------------

global mQuadraturePoly aProdPrimeQuadraturePoly

% Create individual polynomials
mProdQuadraturePoly 		= computeChebyshev(nProd,mQuadratureGridZeros(:,1));
mCapitalQuadraturePoly 	= computeChebyshev(nCapital,mQuadratureGridZeros(:,2));

% Compute tensor product of polynomials
mQuadraturePoly 			= zeros(nStateQuadrature,nState);
for iState = 1:nStateQuadrature
	[vProd,vCapital] 					= ndgrid(mProdQuadraturePoly(iState,:),mCapitalQuadraturePoly(iState,:));
	mQuadraturePoly(iState,:) 	= reshape(vProd .* vCapital,nState,1);
end
clear vProd vCapital

% Compute polynomials over future shocks
aProdPrimeQuadraturePoly = reshape(computeChebyshev(nProd,reshape(mProdPrimeQuadratureZeros,nShocks * ...
										nStateQuadrature,1)),nShocks,nStateQuadrature,nProd);

	
%---------------------------------------------------------------
% Derivative of Chebyshev polynomials over collocation nodes
%		(used to compute first-order condition in dynamic model)	
%---------------------------------------------------------------

global mStatePoly2

% Create individual polynomials
mProdPoly 			= computeChebyshev(nProd,mStateGridZeros(:,1));
mCapitalPoly2 	= computeChebyshev2(nCapital-1,mStateGridZeros(:,2));

% Create tensor product in polynomials
mStatePoly2 		= zeros(nState,nState-nProd);
for iState = 1:nState
	[vProd,vCapital] 			= ndgrid(mProdPoly(iState,:),mCapitalPoly2(iState,:));
	mStatePoly2(iState,:) 	= reshape(vProd .* vCapital,nState-nProd,1);
end
clear mCapitalPoly2 vProd vCapital


%---------------------------------------------------------------
% Derivative of Chebyshev polynomials over fine grid
%		(used to plot marginal value function)	
%---------------------------------------------------------------

global mFinePoly2

% Create individual polynomials
mProdPoly 			= computeChebyshev(nProd,mFineGridZeros(:,1));
mCapitalPoly2 	= computeChebyshev2(nCapital-1,mFineGridZeros(:,2));

% Create tensor product in polynomials
mFinePoly2 		= zeros(nStateFine,nState-nProd);
for iState = 1:nStateFine
	[vProd,vCapital]		 = ndgrid(mProdPoly(iState,:),mCapitalPoly2(iState,:));
	mFinePoly2(iState,:) = reshape(vProd .* vCapital,nState-nProd,1);
end
clear mCapitalPoly2 vProd vCapital
