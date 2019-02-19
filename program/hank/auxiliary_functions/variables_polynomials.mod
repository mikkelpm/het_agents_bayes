// Declare variables for firstOrderDynamics.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Conditional expectation coefficients
//----------------------------------------------------------------

@#for iShare in 1 : nShare
	@#for iPower in 1 : nAssets
		@#for iz in 1 : nz
			var expectationCoefficient_@{iz}_@{iPower}_@{iShare};
		@#endfor
	@#endfor
@#endfor

//----------------------------------------------------------------
// Density of households away from borrowing constraint
//----------------------------------------------------------------

@#for iShare in 1 : nShare
	@#for iMoment in 1 : nMeasure
		@#for iz in 1 : nz
			var moment_@{iz}_@{iMoment}_@{iShare};
		@#endfor
	@#endfor
@#endfor

// Parameters of the distribution
@#for iShare in 1 : nShare
	@#for iParameter in 1 : nMeasure
		@#for iz in 1 : nz
			var measureCoefficient_@{iz}_@{iParameter}_@{iShare};
		@#endfor
	@#endfor
@#endfor

//----------------------------------------------------------------
// Mass at borrowing constraint
//----------------------------------------------------------------

@#for iShare in 1 : nShare
	@#for iz in 1 : nz
		var mHat_@{iz}_@{iShare};
	@#endfor
@#endfor

//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var r i w ppi;

//----------------------------------------------------------------
// Driving variables
//----------------------------------------------------------------

var A uM;

//----------------------------------------------------------------
// Aggregate variables
//----------------------------------------------------------------

var N T d Y C;

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

varexo epsilonA epsilonM;