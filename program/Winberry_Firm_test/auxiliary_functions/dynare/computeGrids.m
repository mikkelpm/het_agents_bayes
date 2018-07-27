% Computes grids used in various approximations
%
% Thomas Winberry, February 14th, 2018

%---------------------------------------------------------------
% Grids to integrate idisoyncratic shocks (innovations to firm-level TFP)
%---------------------------------------------------------------

global vShocksGrid vShocksWeights

[vShocksGrid,vShocksWeights] 	= computeGaussHermiteQuadrature(nShocks);
vShocksGrid 							= sqrt(2) * vShocksGrid;
vShocksWeights 						= (pi ^ (-.5)) * vShocksWeights;


%---------------------------------------------------------------
% Grids for approximating value function and capital accumulation decision
%		conditional on adjusting (Chebyshev collocation nodes)
%---------------------------------------------------------------

global vProdZeros vCapitalZeros mStateGridZeros mStateGrid mProdPrimeZeros

% Compute zeros Chebyshev polynomials
vProdZeros 		= -cos(((2*(1:nProd)-1)' * pi) / (2 * nProd));
vCapitalZeros 	= -cos(((2*(1:nCapital)-1)' * pi) / (2 * nCapital));

% Make tensor product grid
[vProd,vCapital] 	= ndgrid(vProdZeros,vCapitalZeros);
mStateGridZeros = [reshape(vProd,nState,1) reshape(vCapital,nState,1)];

% Scale up to state space
mStateGrid 		= [scaleUp(mStateGridZeros(:,1),prodMin,prodMax) scaleUp(mStateGridZeros(:,2),capitalMin,...
							capitalMax)];

% Make grid of future productivity shocks (for computing expectations)
mProdPrimeZeros = rrhoProd * repmat(mStateGrid(:,1)',[nShocks 1]) + ssigmaProd * repmat(vShocksGrid,[1 nState]);
mProdPrimeZeros = scaleDown(mProdPrimeZeros,prodMin,prodMax);


%---------------------------------------------------------------
% Grids for approximating histogram and plotting functions 
%		(finer than Chebyshev collocation nodes)
%---------------------------------------------------------------

global vProdGridFine vCapitalGridFine mFineGrid mFineGridZeros mProdPrimeFineZeros mProdTransition

% Individual grids
vProdGridFine 		= linspace(prodMin,prodMax,nProdFine)';
vCapitalGridFine 		= linspace(capitalMin,capitalMax,nCapitalFine)';

% Make tensor product grid
[vProd,vCapital] 		= ndgrid(vProdGridFine,vCapitalGridFine);
mFineGrid 				= [reshape(vProd,nStateFine,1) reshape(vCapital,nStateFine,1)];
clear vProd vCapital		% memory management

% Scale down to [-1,1] for polynomials
mFineGridZeros 		= [scaleDown(mFineGrid(:,1),prodMin,prodMax) scaleDown(mFineGrid(:,2),capitalMin,capitalMax)];

% Make grid of future productivity shocks
mProdPrimeFineZeros = rrhoProd * repmat(mFineGrid(:,1)',[nShocks 1]) + ssigmaProd * repmat(vShocksGrid,[1 nStateFine]);
mProdPrimeFineZeros = scaleDown(mProdPrimeFineZeros,prodMin,prodMax);

% Make Tauchen transition matrix for productivity shocks
vPlus 						= zeros(1,nProdFine);
vPlus(1,1:nProdFine-1) = vProdGridFine(2:nProdFine)';
vPlus(1,nProdFine) 		= 1e9;

vMinus 						= zeros(1,nProdFine);
vMinus(1,1) 				= -1e9;
vMinus(1,2:nProdFine) 	= vProdGridFine(1:nProdFine-1)';

mPlusCutoff 				= .5 * ones(nProdFine,1) * (vProdGridFine' + vPlus) - rrhoProd * vProdGridFine * ones(1,nProdFine);
mMinusCutoff 			= .5 * ones(nProdFine,1) * (vProdGridFine' + vMinus) - rrhoProd * vProdGridFine * ones(1,nProdFine);
mProdTransition 			= normcdf(mPlusCutoff,0,ssigmaProd) - normcdf(mMinusCutoff,0,ssigmaProd);


%---------------------------------------------------------------
% Nodes and weights for integrating parametric family
%---------------------------------------------------------------

global vProdQuadrature vCapitalQuadrature mQuadratureGrid vQuadratureWeights mQuadratureGridZeros mProdPrimeQuadrature

% Compute grids in the interval [-1,1]
[vProd,vWeightsProd] 			= computeGaussLegendreQuadrature(nProdQuadrature);
[vCapital,vWeightsCapital] 	= computeGaussLegendreQuadrature(nCapitalQuadrature);

% Scale up grid correctly
vProdQuadrature 				= scaleUp(vProd,prodMin,prodMax);
vCapitalQuadrature 				= scaleUp(vCapital,capitalMin,capitalMax);

% Make tensor product grid
[vProd,vCapital] 					= ndgrid(vProdQuadrature,vCapitalQuadrature);
mQuadratureGrid 				= [reshape(vProd,nStateQuadrature,1) reshape(vCapital,nStateQuadrature,1)];

% Scale down to [-1,1]
mQuadratureGridZeros 		= [scaleDown(mQuadratureGrid(:,1),prodMin,prodMax) scaleDown(mQuadratureGrid(:,2),capitalMin,capitalMax)];

% Make tensor product weights
vWeightsProd 					= ((prodMax - prodMin) / 2) * vWeightsProd;
vWeightsCapital 					= ((capitalMax - capitalMin) / 2) * vWeightsCapital;
[vWeightsProd,vWeightsCapital] = ndgrid(vWeightsProd,vWeightsCapital);
vQuadratureWeights 			= reshape(vWeightsProd .* vWeightsCapital,nStateQuadrature,1);

% Make grid over future productivity shocks (useful in computing decisions)
mProdPrimeQuadrature 		= rrhoProd * repmat(mQuadratureGrid(:,1)',[nShocks 1]) + ssigmaProd * repmat(vShocksGrid,[1 nStateQuadrature]);
mProdPrimeQuadratureZeros = scaleDown(mProdPrimeQuadrature,prodMin,prodMax);
