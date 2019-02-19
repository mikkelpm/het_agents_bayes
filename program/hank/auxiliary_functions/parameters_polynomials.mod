// Declare parameters and load in their values for firstOrderDynamics_polynomials.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Preliminaries
//----------------------------------------------------------------

// Load in files containing parameters
economicParameters = load('economicParameters');
approximationParameters = load('approximationParameters');
grids = load('grids');
polynomials = load('polynomials');

// Define economic parameters 
parameters bbeta ppsi nnu bbBar eepsilon ttau vvarthetaB vvarthetaT
	A_SS w_SS r_RepSS N_RepSS rrhoTFP ssigmaTFP rrhoMP ssigmaMP ttheta pphi;
//Load in their values
@#define nEconomicParameters = 18
for iParam = 1 : @{nEconomicParameters}
	parameterName = deblank(M_.param_names(iParam,:));
	if isfield(economicParameters,parameterName)
		M_.params(iParam) = eval(['economicParameters.' parameterName]);
	end
end

@#define nCounter = nEconomicParameters

// z transition matrix
// Have to impose parameters used as counters directly
@#define nz = 3
@#for iz in 1 : nz
	@#for izPrime in 1 : nz
		parameters zTransition_@{iz}_@{izPrime};
	@#endfor
@#endfor
for iz = 1 : nz
	for izPrime = 1 : nz
		M_.params(@{nCounter} + @{nz} * (iz - 1) + izPrime) = ...
			economicParameters.mzTransition(iz,izPrime);
	end
end

@#define nCounter = nCounter + nz*nz

// Mass of invariant distribution of idiosyncratic shocks 
// Define
@#for iz in 1 : nz
	parameters zMass_@{iz};
@#endfor
// Assign values
for iz = 1 : @{nz}
	M_.params(@{nCounter} + iz) = economicParameters.vzInvariant(iz);
end
@#define nCounter = nCounter + nz

// Equity shares
// Have to impose parameters used as counters directly
@#define nShare = 3
// Profit shares for each household type
// Define
@#for iShare in 1 : nShare
	parameters shareGrid_@{iShare};
@#endfor
// Assign values
for iShare = 1 : @{nShare}
	M_.params(@{nCounter} + iShare) = economicParameters.vShareGrid(iShare);
end
@#define nCounter = nCounter + nShare

// Fractions of each household type
// Define
@#for iShare in 1 : nShare
	parameters shareMass_@{iShare};
@#endfor
// Assign values
for iShare = 1 : @{nShare}
	M_.params(@{nCounter} + iShare) = economicParameters.vShareFraction(iShare);
end
@#define nCounter = nCounter + nShare

// Define some of the approximation parameters
parameters nz nAssets nState nAssetsFine nStateFine
    nAssetsQuadrature nStateQuadrature nMeasure nMeasureCoefficients nShare
    assetsMin assetsMax
    tolerance_SS_root tolerance_SS_invhist maxIterations tolerance dampening numNewton;
// Load in their values
@#define nApproximationParameters = 18
for iParam = 1 : @{nApproximationParameters}
	parameterName = deblank(M_.param_names(@{nCounter} + iParam,:));
	if isfield(approximationParameters,parameterName)
		M_.params(@{nCounter} + iParam) = eval(['approximationParameters.' parameterName]);
	end
end

@#define nCounter = nCounter + nApproximationParameters

// Have to impose parameters used as counters directly
@#define nAssets = 25
@#define nMeasure = 3
@#define nAssetsQuadrature = 8

@#define nState = nz * nAssets
@#define nStateQuadrature = nz * nAssetsQuadrature

//
//----------------------------------------------------------------
// Grids for approximating conditional expectation
//----------------------------------------------------------------

// 
// Idiosyncratic productivity
//

// Define the grids
@#for iz in 1 : nz
	parameters zGrid_@{iz};
@#endfor

// Assign values
for iz = 1 : @{nz}
	M_.params(@{nCounter} + iz) = economicParameters.vzGrid(iz);
end

@#define nCounter = nCounter + nz

//
// Assets
//

// Define the grids
@#for iAssets in 1 : nAssets
	parameters assetsGrid_@{iAssets};
@#endfor

// Assign values (must be in the same order that the parameters were declared)
for iAssets = 1 : @{nAssets}
	M_.params(@{nCounter} + iAssets) = grids.vAssetsGrid(iAssets);
end

// Update counter for future parameter assignments
@#define nCounter = nCounter + nAssets

//----------------------------------------------------------------
// Quadrature grid and weights
//----------------------------------------------------------------

// Define the parameters
@#for iAssets in 1 : nAssetsQuadrature
	parameters quadratureGrid_@{iAssets};
	parameters quadratureWeights_@{iAssets};
@#endfor

// Assign values
for iAssets = 1 : @{nAssetsQuadrature}
	M_.params(@{nCounter} + 2 * (iAssets - 1) + 1) = grids.vAssetsGridQuadrature(iAssets);
	M_.params(@{nCounter} + 2 * (iAssets - 1) + 2) = grids.vQuadratureWeights(iAssets);
end

@#define nCounter = nCounter + 2 * nAssetsQuadrature

//----------------------------------------------------------------
// Conditional expectation polynomials
//----------------------------------------------------------------

//
// Chebyshev polynomials
//

// Define the parameters
@#for iAssets in 1 : nAssets
	@#for iPower in 1 : nAssets
		parameters expectationPoly_@{iAssets}_@{iPower};
	@#endfor
@#endfor

// Assign values
for iAssets = 1 : @{nAssets}
	for iPower = 1 : @{nAssets}
		M_.params(@{nCounter} + @{nAssets} * (iAssets - 1) + iPower) = polynomials.vAssetsPoly(iAssets,iPower);
	end
end

@#define nCounter = nCounter + nAssets * nAssets

// 
// Squared terms in Chebyshev interpolation
//

// Define the parameters
@#for iAssets in 1 : nAssets
	parameters expectationPolySquared_@{iAssets};
@#endfor

// Assign the values
for iAssets = 1 : @{nAssets}
	M_.params(@{nCounter} + iAssets) = polynomials.vAssetsPolySquared(iAssets);	
end

@#define nCounter = nCounter + nAssets

//----------------------------------------------------------------
// Quadrature grid polynomials
//----------------------------------------------------------------

// Define the parameters
@#for iAssets in 1 : nAssetsQuadrature
	@#for iPower in 1 : nAssets
		parameters quadraturePoly_@{iAssets}_@{iPower};		
	@#endfor
@#endfor

// Assign values
for iAssets = 1 : @{nAssetsQuadrature}
	for iPower = 1 : @{nAssets}	
		M_.params(@{nCounter} + @{nAssets} * (iAssets - 1) + iPower) = polynomials.vAssetsPolyQuadrature(iAssets,iPower);		
	end	
end

@#define nCounter = nCounter + nAssetsQuadrature * nAssets

//----------------------------------------------------------------
// Borrowing constraint polynomials
//----------------------------------------------------------------

// Define the parameters
@#for iPower in 1 : nAssets
	parameters bcPoly_@{iPower};
@#endfor

// Assign values
for iPower = 1 : @{nAssets}
	M_.params(@{nCounter} + iPower) = polynomials.vAssetsPolyBC(iPower);
end

@#define nCounter = nCounter + nAssets

//----------------------------------------------------------------
// Parameters depending on steady state (values changed in the steady state file)
//----------------------------------------------------------------

parameters r_SS B_SS G_SS;
//Set to 0 for now to satisfy Dynare
@#define nSSParameters = 3
for iParam = 1 : @{nSSParameters}
    M_.params(@{nCounter} + iParam) = 0;
end