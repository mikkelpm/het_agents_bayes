
%% Load in files containing parameters

economicParameters = load('economicParameters');
approximationParameters = load('approximationParameters');
grids = load('grids');
polynomials = load('polynomials');

%% Economic parameters
nEconomicParameters = 18;
for iParam = 1 : nEconomicParameters
	parameterName = deblank(M_.param_names(iParam,:));
	if isfield(economicParameters,parameterName)
		M_.params(iParam) = eval(['economicParameters.' parameterName]);
	end
end
	
%% Some of the approximation parameters
nApproximationParameters = 13;
for iParam = 1 : nApproximationParameters
	parameterName = deblank(M_.param_names(nEconomicParameters + iParam,:));	
	if isfield(approximationParameters,parameterName)
		M_.params(nEconomicParameters + iParam) = eval(['approximationParameters.' parameterName]);
	end	
end

nCounter = nEconomicParameters + nApproximationParameters;

%% Value function grids

for iState = 1 : nState
	for iVariable = 1 : 2
		M_.params(nCounter + 2 * (iState - 1) + iVariable) = grids.mStateGrid(iState,iVariable);
	end
end

nCounter = nCounter + 2 * nState;

%% Idiosyncratic shocks

% Quadrature grid
for iShock = 1 : nShocks
	M_.params(nCounter + iShock) = grids.vShocksGrid(iShock);
end

nCounter = nCounter + nShocks;

% Quadrature weights
for iShock = 1 : nShocks
	M_.params(nCounter + iShock) = grids.vShocksWeights(iShock);	
end

nCounter = nCounter + nShocks;

%% Grids for integrating measure

% Grid
for iState = 1 : nStateQuadrature
	for iVariable = 1 : 2
		M_.params(nCounter + 2 * (iState - 1) + iVariable) = grids.mQuadratureGrid(iState,iVariable);	
	end
end

nCounter = nCounter + 2 * nStateQuadrature;

% Quadrature weights
for iState = 1 : nStateQuadrature
	M_.params(nCounter + iState) = grids.vQuadratureWeights(iState);
end

nCounter = nCounter + nStateQuadrature;

% Grid over productivity prime
for iState = 1 : nStateQuadrature
	for iShock = 1 : nShocks	
		M_.params(nCounter + nShocks * (iState - 1) + iShock) = grids.mProdPrimeQuadrature(iShock,iState);		
	end	
end
	
nCounter = nCounter + nShocks * nStateQuadrature;

%% Value function polynomials

% Polynomials over state grid
for iState = 1 : nState
	for iPower = 1 : nState	
		M_.params(nCounter + nState * (iState - 1) + iPower) = polynomials.mStatePoly(iState,iPower);		
	end	
end

nCounter = nCounter + nState * nState;

% Polynomials over productivity prime
for iState = 1 : nState	
	for iShock = 1 : nShocks	
		for iPower = 1 : nProd		
			M_.params(nCounter + nProd * nShocks * (iState - 1) + nProd * (iShock - 1) + iPower) ...
                = polynomials.aProdPrimePoly(iShock,iState,iPower);			
		end		
	end	
end

nCounter = nCounter + nState * nShocks * nProd;

% Derivative of value function
for iState = 1 : nState
	for iPower = 1 : nState - nProd
		M_.params(nCounter + (nState - nProd) * (iState - 1) + iPower) ...
            = polynomials.mStatePoly2(iState,iPower);		
	end	
end

nCounter = nCounter + nState * (nState - nProd);
 
% Squared terms of Chebyshev interpolation
for iState = 1 : nState
	M_.params(nCounter + iState) = polynomials.vStatePolySquared(iState);	
end

nCounter = nCounter + nState;

%% Quadrature grid polynomials

% Polynomials over quadrature grid
for iState = 1 : nStateQuadrature
	for iPower = 1 : nState	
		M_.params(nCounter + nState * (iState - 1) + iPower) ...
            = polynomials.mQuadraturePoly(iState,iPower);		
	end	
end

nCounter = nCounter + nStateQuadrature * nState;


% Polynomials over productivity prime
for iState = 1 : nStateQuadrature
	for iShock = 1 : nShocks
		for iPower = 1 : nProd	
			M_.params(nCounter + nProd * nShocks * (iState - 1) + nProd * (iShock - 1) + iPower)...
                = polynomials.aProdPrimeQuadraturePoly(iShock,iState,iPower);			
		end		
	end	
end

nCounter = nCounter + nStateQuadrature * nShocks * nProd;


% Sample moment var-cov matrix
for iCov = 1 : 9
    M_.params(nCounter + iCov) = 0; % Placeholder (will be set in steady state file)
end

nCounter = nCounter + 9;

