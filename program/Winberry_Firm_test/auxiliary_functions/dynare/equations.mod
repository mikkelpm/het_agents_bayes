// Specifies equations for "dynamicModel.mod"
//
// Thomas Winberry, February 15th, 2018


//----------------------------------------------------------------
// Bellman equation (#equations = nState)
//----------------------------------------------------------------

// Preliminary: expand capitalAdjust (which is defined only for productivity) along entire grid (productivity and capital)
@#for iCapital in 1 : nCapital
	@#for iProd in 1 : nProd
		@#define iState = nProd * (iCapital - 1) + iProd
		# capitalAdjustExpanded_@{iState} = capitalAdjust_@{iProd};
	@#endfor
@#endfor


// Define the Bellman equation for each point in iState in the individual state space
@#for iState in 1 : nState

	// STEP 1: COMPUTE EXPECTED VALUE FUNCTION NEXT PERIOD
	// Step 1a: compute expected value function next period, conditional on adjusting
	# adjustExpectedValueFunction_@{iState} = 0
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (
		@#for iPowerCapital in 1 : nCapital
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * (iPowerCapital - 1) + iPowerProd
				+ valueCoefficient_@{power}(+1) * valueFunctionPrimePolys_@{iState}_@{iShock}_@{iPowerProd} * 
					cos((@{iPowerCapital} - 1) * acos(min(max(2 * ((capitalAdjustExpanded_@{iState} - capitalMin) / 
					(capitalMax - capitalMin)) - 1,-1),1)))
			@#endfor
		@#endfor
		)
	@#endfor
	;
	
	// Step 1b: compute constrained expected value function
	# capitalConstrained_@{iState} = min(max(capitalAdjustExpanded_@{iState},
		(1 - ddelta + aaLower) * valueGrid_@{iState}_2),(1 - ddelta + aaUpper) * valueGrid_@{iState}_2);
	# constrainedExpectedValueFunction_@{iState} = 0
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (
		@#for iPowerCapital in 1 : nCapital
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * (iPowerCapital - 1) + iPowerProd
				+ valueCoefficient_@{power}(+1) * valueFunctionPrimePolys_@{iState}_@{iShock}_@{iPowerProd} * 
					cos((@{iPowerCapital} - 1) * acos(min(max(2 * ((capitalConstrained_@{iState} - capitalMin) /
					(capitalMax - capitalMin)) - 1,-1),1)))
			@#endfor
		@#endfor
		)
	@#endfor
	;
	
	// Step 1c: compute adjustment threshold
	# cutoff_@{iState} = min(max((1 / (wage * marginalUtility)) * (-(marginalUtility / exp(aggregateQ)) * (capitalAdjustExpanded_@{iState} - 
        capitalConstrained_@{iState}) + bbeta * (adjustExpectedValueFunction_@{iState} - constrainedExpectedValueFunction_@{iState})),
		0),ppsiCapital);
        
		
	// STEP 2: COMPUTE ENTIRE BELLMAN EQUATION
	// Left hand side: current value function
	0
	@#for iCoefficient in 1 : nState
		+ valueCoefficient_@{iCoefficient} * valueFunctionPolys_@{iState}_@{iCoefficient}
	@#endfor
	// Right hand side: flow profits plus expected value function
	= marginalUtility * (nnu ^ (nnu / (1 - nnu)) - nnu ^ (1 / (1 - nnu))) * ((exp(aggregateTFP) * 
		exp(valueGrid_@{iState}_1)) ^ (1 / (1 - nnu))) * (valueGrid_@{iState}_2 ^ (ttheta / (1 - nnu))) * 
		(wage ^ (-nnu / (1 - nnu)))	+ (cutoff_@{iState} / ppsiCapital) * 
        (-(marginalUtility / exp(aggregateQ)) * (capitalAdjustExpanded_@{iState} - (1 - ddelta) * valueGrid_@{iState}_2) - 
        (cutoff_@{iState} / 2) * wage * marginalUtility + bbeta * adjustExpectedValueFunction_@{iState}) + 
        (1 - (cutoff_@{iState} / ppsiCapital)) * (-(marginalUtility / exp(aggregateQ)) * (capitalConstrained_@{iState} -
        (1 - ddelta) * valueGrid_@{iState}_2) + bbeta * constrainedExpectedValueFunction_@{iState});
		
@#endfor


//----------------------------------------------------------------
// FOC for adjust capital decision (#equations = nProd)
//----------------------------------------------------------------		

@#for iProd in 1 : nProd
	(marginalUtility / exp(aggregateQ)) = bbeta * (
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (0
		@#for iPowerCapital in 1 : nCapital - 1
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * iPowerCapital + iPowerProd
					+ (2 / (capitalMax - capitalMin)) * @{iPowerCapital} * valueCoefficient_@{power}(+1) * valueFunctionPrimePolys_@{iProd}_@{iShock}_@{iPowerProd} *
						(sin((@{iPowerCapital}) * acos(min(max(2 * ((capitalAdjust_@{iProd} - capitalMin) / (capitalMax - capitalMin)) - 1, - 1), 1))) / 
						sin(acos(min(max(2 * ((capitalAdjust_@{iProd} - capitalMin) / (capitalMax - capitalMin)) - 1, - 1), 1))))
			@#endfor
		@#endfor
		)
	@#endfor
	);
@#endfor


//----------------------------------------------------------------
// ASIDE: compute objects over quadrature grid for integrating distribution
//		(necessary for the remaining equilibrium conditions, so its convenient 
//		 to do all at once here)
//----------------------------------------------------------------

// Compute coefficients of Chebyshev interpolation of capital accumulation decision, conditional on adjusting
// Will use this interpolation to evaluate capitalAdjust over the entire grid for integrating distribution
@#for iState in 1 : nState
	# capitalAdjustCoefficient_@{iState} = 0
	@#for iCoefficient in 1 : nState
		+ capitalAdjustExpanded_@{iCoefficient} * valueFunctionPolys_@{iCoefficient}_@{iState} / valueFunctionPolySquared_@{iState}
	@#endfor
	;
@#endfor


// Other objects over the grid
@#for iState in 1 : nStateQuadrature

	// PDF of distribution
	# measurePDF_@{iState} = exp(0 + measureCoefficient_1 * (quadratureGrid_@{iState}_1 - moment_1(-1)) +
	measureCoefficient_2 * (log(quadratureGrid_@{iState}_2) - moment_2(-1))
	@#define nIndexCounter = 3
	@#for iOrder in 2 : nMeasure
		@#for iPower in 0 : iOrder
			+ measureCoefficient_@{nIndexCounter} * (((quadratureGrid_@{iState}_1 - moment_1(-1)) ^ (@{iOrder} - @{iPower})) * 
			((log(quadratureGrid_@{iState}_2) - moment_2(-1)) ^ @{iPower}) - moment_@{nIndexCounter}(-1))
			@#define nIndexCounter = nIndexCounter + 1
		@#endfor
	@#endfor
	);
	
	// Capital accumulation decision, conditional on adjusting
	# capitalAdjustQuadrature_@{iState} = min(max(0 + 
	@#for iCoefficient in 1 : nState
		+ capitalAdjustCoefficient_@{iCoefficient} * quadraturePolys_@{iState}_@{iCoefficient}
	@#endfor
	,capitalMin),capitalMax);
	
	// Capital accumulation decision, conditional on not adjusting
	# capitalConstrainedQuadrature_@{iState} = min(max(min(max(capitalAdjustQuadrature_@{iState},(1 - ddelta + aaLower) * quadratureGrid_@{iState}_2),
		(1 - ddelta + aaUpper) * quadratureGrid_@{iState}_2),capitalMin),capitalMax);
		
	// Fixed cost threshold for adjusting or not
	# cutoffQuadrature_@{iState} = min(max((1 / (wage * marginalUtility)) * (-(marginalUtility / exp(aggregateQ)) * (capitalAdjustQuadrature_@{iState} - 
        capitalConstrainedQuadrature_@{iState}) + bbeta * (0
	@#for iShock in 1 : nShocks
		+ shocksWeights_@{iShock} * (
		@#for iPowerCapital in 1 : nCapital
			@#for iPowerProd in 1 : nProd
				@#define power = nProd * (iPowerCapital - 1) + iPowerProd
				+ valueCoefficient_@{power}(+1) * quadraturePrimePolys_@{iState}_@{iShock}_@{iPowerProd} * (
					cos((@{iPowerCapital} - 1) * acos(min(max(2 * ((capitalAdjustQuadrature_@{iState} - capitalMin) 
					/ (capitalMax - capitalMin)) - 1, - 1), 1))) - cos((@{iPowerCapital} - 1) * acos(min(max(2 * 
					((capitalConstrainedQuadrature_@{iState} - capitalMin) / (capitalMax - capitalMin)) - 1, - 1), 1))))
			@#endfor
		@#endfor
		)
	@#endfor
	)),0),ppsiCapital);
	
@#endfor


// Compute total mass of distribution for normalization
# totalMass = 0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState}
@#endfor
;


//----------------------------------------------------------------
// Relationship between moments of distribution and parameters
// (#equations = nMeasureCoefficients)
//----------------------------------------------------------------

// First moments (uncentered)
moment_1(-1) = (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * quadratureGrid_@{iState}_1 * measurePDF_@{iState}
@#endfor
) / totalMass;

moment_2(-1) = (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * log(quadratureGrid_@{iState}_2) * measurePDF_@{iState}
@#endfor
) / totalMass;
	
	
// Higher order moments (centered)
@#define nIndexCounter = 3
@#for iOrder in 2 : nMeasure
	@#for iPower in 0 : iOrder
		moment_@{nIndexCounter}(-1) = (0
		@#for iState in 1 : nStateQuadrature
			+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((quadratureGrid_@{iState}_1 - moment_1(-1)) ^ (@{iOrder} - @{iPower})) * 
			((log(quadratureGrid_@{iState}_2) - moment_2(-1)) ^ @{iPower})
		@#endfor
		) / totalMass;
		@#define nIndexCounter = nIndexCounter + 1
	@#endfor
@#endfor


//----------------------------------------------------------------
// Law of motion for distribution 
// (#equations = nMeasureCoefficients)
//----------------------------------------------------------------

// First moment of productivity (uncentered)
moment_1 = (0
@#for iShock in 1 : nShocks
	+ shocksWeights_@{iShock} * (
	@#for iState in 1 : nStateQuadrature
		+ quadratureWeights_@{iState} * measurePDF_@{iState} * quadratureProdPrimeGrid_@{iState}_@{iShock}
	@#endfor
	)
@#endfor
) / totalMass;

// First moment of capital (uncentered)
moment_2 = (0
@#for iShock in 1 : nShocks
	+ shocksWeights_@{iShock} * (
	@#for iState in 1 : nStateQuadrature
		+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((cutoffQuadrature_@{iState} / ppsiCapital) * log(capitalAdjustQuadrature_@{iState}) + 
			(1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * log(capitalConstrainedQuadrature_@{iState}))
	@#endfor
	)
@#endfor
) / totalMass;

// Higher order moments (centered)
@#define nIndexCounter = 3
@#for iOrder in 2 : nMeasure
	@#for iPower in 0 : iOrder
		moment_@{nIndexCounter} = (0
		@#for iShock in 1 : nShocks
			+ shocksWeights_@{iShock} * (
			@#for iState in 1 : nStateQuadrature
				+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((quadratureProdPrimeGrid_@{iState}_@{iShock} - moment_1) ^ (@{iOrder} - @{iPower})) * 
					(((cutoffQuadrature_@{iState} / ppsiCapital) * (log(capitalAdjustQuadrature_@{iState}) - moment_2) ^ @{iPower}) + 
					((1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * (log(capitalConstrainedQuadrature_@{iState}) - moment_2) ^ @{iPower}))
			@#endfor
			)
		@#endfor
		) / totalMass;
		@#define nIndexCounter = nIndexCounter + 1
	@#endfor
@#endfor

//----------------------------------------------------------------
// Labor Market clearing (# equations = 2)
//----------------------------------------------------------------

// Definition of aggregate hours from labor demand
aggregateHours = (0 
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState} * (((nnu * exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) 
		* (quadratureGrid_@{iState}_2 ^ ttheta)) / wage) ^ (1 / (1 - nnu)) + ((cutoffQuadrature_@{iState} ^ 2) / (2 * ppsiCapital)))
@#endfor
) / totalMass;

// Labor demand = labor supply
aggregateHours = ((wage * marginalUtility) / cchi) ^ (1 / pphi);


//----------------------------------------------------------------
// Output Market clearing (# equations = 2)
//----------------------------------------------------------------

// Definition of aggregate consumption
aggregateConsumption = (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState} * (exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) * (quadratureGrid_@{iState}_2 ^ ttheta) * 
		((((nnu * exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) * (quadratureGrid_@{iState}_2 ^ ttheta)) / wage) ^ 
		(1 / (1 - nnu))) ^ nnu) + (1 / exp(aggregateQ)) * (1 - ddelta) * quadratureGrid_@{iState}_2 + (cutoffQuadrature_@{iState} / ppsiCapital) * 
		(-(1 / exp(aggregateQ)) * capitalAdjustQuadrature_@{iState}) + 
		(1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * (-(1 / exp(aggregateQ)) * capitalConstrainedQuadrature_@{iState}))
@#endfor
) / totalMass;

// Marginal utility = u'(aggregate consumption)
marginalUtility = aggregateConsumption ^ (-ssigma);


//----------------------------------------------------------------
// Law of motion for aggregate shocks (# equations = 2)
//----------------------------------------------------------------

aggregateTFP = rrhoTFP * aggregateTFP(-1) + ssigmaTFP * aggregateTFPShock;
aggregateQ = rrhoQ * aggregateQ(-1) + ssigmaQ * aggregateQShock + corrTFPQ * aggregateTFPShock;


//----------------------------------------------------------------
// Auxiliary variables of interest (# equations = 11)
//----------------------------------------------------------------

// Definition of log aggregate output
logAggregateOutput = log((0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState} * (exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) * (quadratureGrid_@{iState}_2 ^ ttheta) * 
		((((nnu * exp(aggregateTFP) * exp(quadratureGrid_@{iState}_1) * (quadratureGrid_@{iState}_2 ^ ttheta)) / wage) ^ 
		(1 / (1 - nnu))) ^ nnu))	
@#endfor
) / totalMass);

// Definition of log aggregate investment
logAggregateInvestment = log((1 / exp(aggregateQ)) * (0
@#for iState in 1 : nStateQuadrature
	+ quadratureWeights_@{iState} * measurePDF_@{iState} * ((cutoffQuadrature_@{iState} / ppsiCapital) * capitalAdjustQuadrature_@{iState} + 
		(1 - (cutoffQuadrature_@{iState} / ppsiCapital)) * capitalConstrainedQuadrature_@{iState} - (1 - ddelta) * quadratureGrid_@{iState}_2)
@#endfor
) / totalMass);


// Logs of aggregate variables already defined
logAggregateHours = log(aggregateHours);
logAggregateConsumption = log(aggregateConsumption);
logWage = log(wage);
logMarginalUtility	= log(marginalUtility);

// Risk-free interest rate
expectedMarginalUtilityPrime = marginalUtility(+1);
realInterestRate = 100 * ((marginalUtility / (bbeta * expectedMarginalUtilityPrime)) - 1);
