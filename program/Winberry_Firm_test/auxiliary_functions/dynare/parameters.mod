// Load and declare parameters for "dynamicModel.mod"
// NB: parameters include the grids and polynomials defined over the grids
//
// Thomas Winberry, February 15th, 2018


//----------------------------------------------------------------
// Preliminaries
//----------------------------------------------------------------

// Load in files containing parameters
economicParameters 			= load('economicParameters');
approximationParameters 	= load('approximationParameters');
grids 								= load('grids');
polynomials 						= load('polynomials');

// Define economic parameters
parameters ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ cchi;

// Load in their values
@#define nEconomicParameters = 18
for iParam = 1 : @{nEconomicParameters}
	parameterName = deblank(M_.param_names(iParam,:));
	if isfield(economicParameters,parameterName)
		M_.params(iParam) = eval(['economicParameters.' parameterName]);
	end
end
	
// Define some of the approximation parameters
parameters nProd nCapital nState nProdQuadrature nCapitalQuadrature nStateQuadrature nMeasureCoefficients
	nMeasure prodMin prodMax capitalMin capitalMax nShocks;
@#define nApproximationParameters = 13
for iParam = 1 : @{nApproximationParameters}
	parameterName = deblank(M_.param_names(@{nEconomicParameters} + iParam,:));	
	if isfield(approximationParameters,parameterName)
		M_.params(@{nEconomicParameters} + iParam) = eval(['approximationParameters.' parameterName]);
	end	
end

// Have to impose the parameters used as counters below exogenously
@#define nProd 						= 3
@#define nCapital 					= 5
//@#define nMeasure 				= 2
@#define nShocks 					= 3
@#define nProdQuadrature 		= 8
@#define nCapitalQuadrature 	= 10
@#define nState 						= nProd * nCapital
@#define nStateQuadrature 		= nProdQuadrature * nCapitalQuadrature
@#define nMeasureCoefficients 	= (nMeasure * (nMeasure + 1)) / 2 + nMeasure


//
//----------------------------------------------------------------
// Value function grids
//----------------------------------------------------------------

@#define nCounter = nEconomicParameters + nApproximationParameters

// Define the grids
@#for iState in 1 : nState
    @#for iVariable in 1 : 2
        parameters valueGrid_@{iState}_@{iVariable};
    @#endfor
@#endfor

// Assign values
for iState = 1 : @{nState}
	for iVariable = 1 : 2
		M_.params(@{nCounter} + 2 * (iState - 1) + iVariable) = grids.mStateGrid(iState,iVariable);
	end
end

// Index of the parameter vector
@#define nCounter = nCounter + 2 * nState


//----------------------------------------------------------------
// Idiosyncratic shocks
//----------------------------------------------------------------

//
// Quadrature grid
//

// Define the parameters
@#for iShock in 1 : nShocks
	parameters shocksGrid_@{iShock};	
@#endfor
	
// Assign values
for iShock = 1 : @{nShocks}
	M_.params(@{nCounter} + iShock) = grids.vShocksGrid(iShock);
end

@#define nCounter = nCounter + nShocks

//
// Quadrature weights
//

// Define the parameters
@#for iShock in 1 : nShocks  
    parameters shocksWeights_@{iShock};  
@#endfor

// Assign values
for iShock = 1 : @{nShocks}
	M_.params(@{nCounter} + iShock) = grids.vShocksWeights(iShock);	
end

@#define nCounter = nCounter + nShocks


//----------------------------------------------------------------
// Grids for integrating measure
//----------------------------------------------------------------

// 
// Grid
//

// Define the parameters
@#for iState in 1 : nStateQuadrature
	@#for iVariable in 1 : 2	
		parameters quadratureGrid_@{iState}_@{iVariable};	
	@#endfor
@#endfor

// Assign values
for iState = 1 : @{nStateQuadrature}
	for iVariable = 1 : 2
		M_.params(@{nCounter} + 2 * (iState - 1) + iVariable) = grids.mQuadratureGrid(iState,iVariable);	
	end
end

@#define nCounter = nCounter + 2 * nStateQuadrature


//
// Quadrature weights
//

// Define the parameters
@#for iState in 1 : nStateQuadrature
	parameters quadratureWeights_@{iState};	
@#endfor

// Assign values
for iState = 1 : @{nStateQuadrature}
	M_.params(@{nCounter} + iState) = grids.vQuadratureWeights(iState);
end

@#define nCounter = nCounter + nStateQuadrature


//
// Grid over productivity prime
//

// Define the parameters
@#for iState in 1 : nStateQuadrature
	@#for iShock in 1 : nShocks
		parameters quadratureProdPrimeGrid_@{iState}_@{iShock};
	@#endfor	
@#endfor

// Assign values
for iState = 1 : @{nStateQuadrature}
	for iShock = 1 : @{nShocks}	
		M_.params(@{nCounter} + @{nShocks} * (iState - 1) + iShock) = grids.mProdPrimeQuadrature(iShock,iState);		
	end	
end
	
@#define nCounter = nCounter + nShocks * nStateQuadrature	


//----------------------------------------------------------------
// Value function polynomials
//----------------------------------------------------------------

// 
// Polynomials over state grid
//

// Define the parameters
@#for iState in 1 : nState
	@#for iPower in 1 : nState	
		parameters valueFunctionPolys_@{iState}_@{iPower};		
	@#endfor
@#endfor

// Assign values
for iState = 1 : @{nState}
	for iPower = 1 : @{nState}	
		M_.params(@{nCounter} + @{nState} * (iState - 1) + iPower) = polynomials.mStatePoly(iState,iPower);		
	end	
end

@#define nCounter = nCounter + nState * nState


//
// Polynomials over productivity prime
//

// Define the parameters
@#for iState in 1 : nState	
	@#for iShock in 1 : nShocks	
		@#for iPower in 1 : nProd			
			parameters valueFunctionPrimePolys_@{iState}_@{iShock}_@{iPower};		
		@#endfor	
	@#endfor	
@#endfor

// Assign values
for iState = 1 : @{nState}	
	for iShock = 1 : @{nShocks}	
		for iPower = 1 : @{nProd}		
			M_.params(@{nCounter} + @{nProd} * @{nShocks} * (iState - 1) + @{nProd} * (iShock - 1) + iPower) = polynomials.aProdPrimePoly(iShock,iState,iPower);			
		end		
	end	
end

@#define nCounter = nCounter + nState * nShocks * nProd


// 
// Derivative of value function
//

// Define the parameters
@#for iState in 1 : nState
	@#for iPower in 1 : nState - nProd	
		parameters marginalValueFunctionPolys_@{iState}_@{iPower};		
	@#endfor
@#endfor

// Assign values
for iState = 1 : @{nState}
	for iPower = 1 : @{nState} - @{nProd}
		M_.params(@{nCounter} + (@{nState} - @{nProd}) * (iState - 1) + iPower) = polynomials.mStatePoly2(iState,iPower);		
	end	
end

@#define nCounter = nCounter + nState * (nState - nProd)


// 
// Squared terms of Chebyshev interpolation
//

// Define the parameters
@#for iState in 1 : nState
	parameters valueFunctionPolySquared_@{iState};
@#endfor

// Assign the values
for iState = 1 : @{nState}
	M_.params(@{nCounter} + iState) = polynomials.vStatePolySquared(iState);	
end

@#define nCounter = nCounter + nState


//----------------------------------------------------------------
// Quadrature grid polynomials
//----------------------------------------------------------------

// 
// Polynomials over quadrature grid
//

// Define the parameters
@#for iState in 1 : nStateQuadrature
	@#for iPower in 1 : nState	
		parameters quadraturePolys_@{iState}_@{iPower};		
	@#endfor
@#endfor

// Assign values
for iState = 1 : @{nStateQuadrature}
	for iPower = 1 : @{nState}	
		M_.params(@{nCounter} + @{nState} * (iState - 1) + iPower) = polynomials.mQuadraturePoly(iState,iPower);		
	end	
end

@#define nCounter = nCounter + nStateQuadrature * nState


//
// Polynomials over productivity prime
//

// Define the parameters
@#for iState in 1 : nStateQuadrature	
	@#for iShock in 1 : nShocks
		@#for iPower in 1 : nProd			
			parameters quadraturePrimePolys_@{iState}_@{iShock}_@{iPower};	
		@#endfor	
	@#endfor	
@#endfor

// Assign values
for iState = 1 : @{nStateQuadrature}
	for iShock = 1 : @{nShocks}
		for iPower = 1 : @{nProd}		
			M_.params(@{nCounter} + @{nProd} * @{nShocks} * (iState - 1) + @{nProd} * (iShock - 1) + iPower) = polynomials.aProdPrimeQuadraturePoly(iShock,iState,iPower);			
		end		
	end	
end

@#define nCounter = nCounter + nStateQuadrature * nShocks * nProd
