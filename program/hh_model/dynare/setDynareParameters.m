
%% Load in files containing parameters

economicParameters = load_mat('economicParameters');
approximationParameters = load_mat('approximationParameters');
grids = load_mat('grids');
polynomials = load_mat('polynomials');


%% Economic parameters

nEconomicParameters = 11; % Note: account for mu_l

for iParam = 1 : nEconomicParameters
	parameterName = deblank(M_.param_names{iParam});
	if isfield(economicParameters,parameterName)
		M_.params(iParam) = economicParameters.(parameterName);
	end
end


%% Epsilon transition matrix

for iEpsilon = 1 : 2
	for iEpsilonPrime = 1 : 2
		M_.params(nEconomicParameters + 2 * (iEpsilon - 1) + iEpsilonPrime) = ...
			economicParameters.mEpsilonTransition(iEpsilon,iEpsilonPrime);
	end
end


%% Mass of invariant distrbution of idiosyncratic shocks

set_param_value('epsilonMass_1', 1 - aggEmployment);
set_param_value('epsilonMass_2', aggEmployment);


%% Some of the approximation parameters 

nApproximationParameters = 15;
for iParam = 1 : nApproximationParameters
	parameterName = deblank(M_.param_names{nEconomicParameters + 6 + iParam});
	if isfield(approximationParameters,parameterName)
		M_.params(nEconomicParameters + 6 + iParam) = approximationParameters.(parameterName);
	end
end

nCounter = nEconomicParameters + 6 + nApproximationParameters;


%% Grids for approximating conditional expectation
 
% Employment

set_param_value('epsilonGrid_1', 0);
set_param_value('epsilonGrid_2', 1);

nCounter = nCounter + 2;

%% Assets

for iAssets = 1 : nAssets
	M_.params(nCounter + iAssets) = grids.vAssetsGrid(iAssets);
end

% Update counter for future parameter assignments
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


%% Quadrature grid polynomials
for iAssets = 1 : nAssetsQuadrature
	for iPower = 1 : nAssets
		M_.params(nCounter + nAssets * (iAssets - 1) + iPower) = polynomials.vAssetsPolyQuadrature(iAssets,iPower);		
	end	
end

nCounter = nCounter + nAssetsQuadrature * nAssets;


%% Borrowing constraint polynomials
for iPower = 1 : nAssets
	M_.params(nCounter + iPower) = polynomials.vAssetsPolyBC(iPower);
end
