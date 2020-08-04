% Sets parameter values 
%
% Thomas Winberry, July 26th, 2016


global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP ssigmaMeas mu_l;

global nEpsilon nAssets nAssetsFine nAssetsQuadrature ...
	nMeasure maxIterations tolerance dampening splineOpt displayOpt;


%% Economic parameters

global mEpsilonTransition vEpsilonInvariant ttau;

% Idioynscratic Shocks

mEpsilonTransition = [uDuration / (1 + uDuration), 1 - (uDuration / (1 + uDuration));...
					 ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration))),...
					 1 - ((1 - aggEmployment) / aggEmployment) * (1 - (uDuration / (1 + uDuration)))];
vEpsilonInvariant = [1 - aggEmployment;aggEmployment];

% Unemployment benefits
ttau = mmu * (1 - aggEmployment) / aggEmployment;


%% Approximation parameters

global nState assetsMin assetsMax nStateFine nStateQuadrature nMeasureCoefficients kRepSS;

% Approximation parameters
nState = nEpsilon * nAssets;

% Bounds on grid space
kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - aalpha));
assetsMin = aaBar;	assetsMax = 3 * kRepSS;

% Finer grid for analyzing policy functions
nStateFine = nEpsilon * nAssetsFine;

% Approximation of distribution
nStateQuadrature = nEpsilon * nAssetsQuadrature;
nMeasureCoefficients = nEpsilon * nMeasure;


%% Compute approximation tools

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end
