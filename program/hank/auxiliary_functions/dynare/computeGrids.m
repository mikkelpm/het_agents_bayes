% Creates grids to use in various approximations
% 
% Thomas Winberry, July 26th, 2016

%---------------------------------------------------------------
% Grids for approximating individual decisions
%---------------------------------------------------------------


% Zeros of chebyshev polynomial
vAssetsGridZeros = -cos(((2 * (1:nAssets)-1)' * pi) / (2 * nAssets));

% Scale up to state space
vAssetsGrid = scaleUp(vAssetsGridZeros,assetsMin,assetsMax);

	
% Make matrix versions of the grids
mzGrid = repmat(vzGrid,[1 nAssets]);
mAssetsGrid = repmat(vAssetsGrid',[nz 1]);
mzPrimeGrid = repmat(vzGrid,[1 nz*nAssets]);

%---------------------------------------------------------------
% Fine grid, for histogram and plotting functions
%---------------------------------------------------------------

% global vAssetsGridFine vAssetsGridFineZeros mzGridFine mAssetsGridFine mzPrimeGridFine 

% Assets grid
vAssetsGridFine = linspace(assetsMin,assetsMax,nAssetsFine)';

% Scale down to [-1,1]
vAssetsGridFineZeros = scaleDown(vAssetsGridFine,assetsMin,assetsMax);

% Make matrix versions of grids
mzGridFine = repmat(vzGrid,[1 nAssetsFine]);
mAssetsGridFine = repmat(vAssetsGridFine',[nz 1]);
mzPrimeGridFine = repmat(vzGrid,[1 nStateFine]);

%---------------------------------------------------------------
% Quadrature grid, to integrate density (away from borrowing constraint)
%---------------------------------------------------------------

% Compute grid in the interval [-1, 1]
[vAssetsGridQuadratureZeros,vQuadratureWeights] = computeGaussLegendreQuadrature(nAssetsQuadrature);

% Scale up grid
vAssetsGridQuadrature = scaleUp(vAssetsGridQuadratureZeros,assetsMin+1e-1,assetsMax);

% Make matrix versions of the grids
mzGridQuadrature = repmat(vzGrid,[1 nAssetsQuadrature]);
mAssetsGridQuadrature = repmat(vAssetsGridQuadrature', [nz 1]);
