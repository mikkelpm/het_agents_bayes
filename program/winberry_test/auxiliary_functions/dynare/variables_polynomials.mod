// Declare variables for firstOrderDynamics.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Conditional expectation coefficients
//----------------------------------------------------------------

@#for iPower in 1 : nAssets
	var expectationCoefficient_1_@{iPower} expectationCoefficient_2_@{iPower};
@#endfor

//----------------------------------------------------------------
// Density of households away from borrowing constraint
//----------------------------------------------------------------

// Moments of the distribution
@#for iEpsilon in 1 : nEpsilon
    var lag_mHat_@{iEpsilon};
    @#for iMoment in 1 : nMeasure
        var moment_@{iEpsilon}_@{iMoment} smpl_m@{iEpsilon}@{iMoment} lag_moment_@{iEpsilon}_@{iMoment};
    @#endfor
@#endfor

// Parameters of the distribution
@#for iParameter in 1 : nMeasure
    var measureCoefficient_1_@{iParameter} measureCoefficient_2_@{iParameter};  
@#endfor

//----------------------------------------------------------------
// Mass at borrowing constraint
//----------------------------------------------------------------

var mHat_1 mHat_2;

//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var r w;

//----------------------------------------------------------------
// Aggregate TFP
//----------------------------------------------------------------

var aggregateTFP;

//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var logAggregateOutput logAggregateInvestment logAggregateConsumption logWage;

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

varexo aggregateTFPShock;