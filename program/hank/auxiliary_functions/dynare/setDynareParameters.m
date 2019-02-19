
%% Load in files containing parameters

economicParameters = load('economicParameters');
approximationParameters = load('approximationParameters');


%% Economic parameters

nEconomicParameters = 18;
for iParam = 1 : nEconomicParameters
	parameterName = deblank(M_.param_names(iParam,:));
	if isfield(economicParameters,parameterName)
		M_.params(iParam) = eval(['economicParameters.' parameterName]);
	end
end

nCounter = nEconomicParameters;

%% z transition matrix
for iz = 1 : nz
	for izPrime = 1 : nz
		M_.params(nCounter + nz * (iz - 1) + izPrime) = ...
			economicParameters.mzTransition(iz,izPrime);
	end
end

nCounter = nCounter + nz*nz;

%% Mass of invariant distribution of idiosyncratic shocks 
for iz = 1 : nz
	M_.params(nCounter + iz) = economicParameters.vzInvariant(iz);
end

nCounter = nCounter + nz;

%% Equity shares
for iShare = 1 : nShare
	M_.params(nCounter + iShare) = economicParameters.vShareGrid(iShare);
end

nCounter = nCounter + nShare;

%% Fractions of each household type
for iShare = 1 : nShare
	M_.params(nCounter + iShare) = economicParameters.vShareFraction(iShare);
end

nCounter = nCounter + nShare;

%% Define some of the approximation parameters
nApproximationParameters = 18;
for iParam = 1 : nApproximationParameters
	parameterName = deblank(M_.param_names(nCounter + iParam,:));
	if isfield(approximationParameters,parameterName)
		M_.params(nCounter + iParam) = eval(['approximationParameters.' parameterName]);
	end
end

nCounter = nCounter + nApproximationParameters;


%% Grids for approximating conditional expectation

% Idiosyncratic productivity
for iz = 1 : nz
	M_.params(nCounter + iz) = economicParameters.vzGrid(iz);
end

nCounter = nCounter + nz;

% Assets
for iAssets = 1 : nAssets
	M_.params(nCounter + iAssets) = grids.vAssetsGrid(iAssets);
end

nCounter = nCounter + nAssets;


% Quadrature grid and weights
for iAssets = 1 : nAssetsQuadrature
	M_.params(nCounter + 2 * (iAssets - 1) + 1) = grids.vAssetsGridQuadrature(iAssets);
	M_.params(nCounter + 2 * (iAssets - 1) + 2) = grids.vQuadratureWeights(iAssets);
end

nCounter = nCounter + 2 * nAssetsQuadrature;


%% Conditional expectation polynomials
% Chebyshev polynomials
for iAssets = 1 : nAssets
	for iPower = 1 : nAssets
		M_.params(nCounter + nAssets * (iAssets - 1) + iPower) = polynomials.vAssetsPoly(iAssets,iPower);
	end
end

nCounter = nCounter + nAssets * nAssets;


% Squared terms in Chebyshev interpolation
for iAssets = 1 : nAssets
	M_.params(nCounter + iAssets) = polynomials.vAssetsPolySquared(iAssets);	
end

nCounter = nCounter + nAssets;


% Quadrature grid polynomials
for iAssets = 1 : nAssetsQuadrature
	for iPower = 1 : nAssets	
		M_.params(nCounter + nAssets * (iAssets - 1) + iPower) = polynomials.vAssetsPolyQuadrature(iAssets,iPower);		
	end	
end

nCounter = nCounter + nAssetsQuadrature * nAssets;


% Borrowing constraint polynomials
for iPower = 1 : nAssets
	M_.params(nCounter + iPower) = polynomials.vAssetsPolyBC(iPower);
end

nCounter = nCounter + nAssets;


%% Parameters depending on steady state (values changed in the steady state file)

nSSParameters = 3;
for iParam = 1 : nSSParameters
    M_.params(nCounter + iParam) = 0;
end
