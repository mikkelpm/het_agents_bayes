% Sets parameter values 
%
% Thomas Winberry,February 14th, 2018

%% Set model parameters

global ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ  cchi
% 
% %% Set parameters governing approximations
% 
% global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
% 	maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
% 	nCapitalQuadrature kRepSS wRepSS
% 
% 
% %% Grids
% global vShocksGrid vShocksWeights...
% vProdZeros vCapitalZeros mStateGridZeros mStateGrid mProdPrimeZeros...
% vProdGridFine vCapitalGridFine mFineGrid mFineGridZeros mProdPrimeFineZeros mProdTransition...
% vProdQuadrature vCapitalQuadrature mQuadratureGrid vQuadratureWeights mQuadratureGridZeros mProdPrimeQuadrature
% 
% %% Polynomials
% global mStatePoly vStatePolySquared aProdPrimePoly...
% mFinePoly aProdPrimeFinePoly...
% mQuadraturePoly aProdPrimeQuadraturePoly...
% mStatePoly2...
% mFinePoly2

%% Set approximation parameters

global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
	maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
	nCapitalQuadrature kRepSS wRepSS

% Bounds on grid space
prodMin 		= -3 * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
prodMax 		= 3 * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
capitalMin 		= .1 * (exp(prodMin) ^ (1 / (1 - ttheta))) * kRepSS;
capitalMax 		= 2.5 * (exp(prodMax) ^ (1 / (1 - ttheta))) * kRepSS;

%% Compute approximation tools

% Grids
computeGrids;

% Polynomials over grids 
computePolynomials;
