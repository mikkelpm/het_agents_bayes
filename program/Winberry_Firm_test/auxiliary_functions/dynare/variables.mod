// Declare variables for "dynamicModel.mod"
//
// Thomas Winberry, Feburary 15th, 2018

//----------------------------------------------------------------
// Value function coefficients
//----------------------------------------------------------------

@#for iState in 1 : nState
    var valueCoefficient_@{iState};
@#endfor


//----------------------------------------------------------------
// Adjust policy function along productivity grid
//----------------------------------------------------------------

@#for iState in 1 : nProd
    var capitalAdjust_@{iState};
@#endfor


//----------------------------------------------------------------
// Moments
//----------------------------------------------------------------

@#for iMoment in 1 : nMeasureCoefficients
    var moment_@{iMoment}; 
@#endfor


//----------------------------------------------------------------
// Distribution parameters
//----------------------------------------------------------------

@#for iParameter in 1 : nMeasureCoefficients
    var measureCoefficient_@{iParameter};  
@#endfor


//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var wage marginalUtility;


//----------------------------------------------------------------
// Aggregate shocks
//----------------------------------------------------------------

var aggregateTFP aggregateQ;


//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var aggregateConsumption aggregateHours expectedMarginalUtilityPrime realInterestRate
	logAggregateOutput logAggregateConsumption logAggregateInvestment logAggregateHours logWage
	logMarginalUtility;

	
//----------------------------------------------------------------
// Innovations to aggregate shocks
//----------------------------------------------------------------

varexo aggregateTFPShock aggregateQShock;
