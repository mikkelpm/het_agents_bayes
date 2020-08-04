// Solve for aggregate dynamics using DYNARE
//
// Thomas Winberry, February 14th, 2018

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters.mod"


//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables.mod"


//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

	@#include "equations.mod"

end;

varobs logAggregateOutput logAggregateInvestment
@#for iMoment in 1 : nMeasureCoefficients
 smpl_m@{iMoment} 
@#endfor
;
//----------------------------------------------------------------
// Set options
//----------------------------------------------------------------

// Specify shock process (=1 to include shock; =0 to not include shock)
shocks;
    var aggregateTFPShock 	= 1;
    var aggregateQShock   	= 1;
    // Placeholder values for measurement error
    var smpl_m1 = 1;
    var smpl_m2 = 1;
    var smpl_m3 = 1;
    var smpl_m4 = 1;
    var smpl_m5 = 1;
    var smpl_m1, smpl_m2 = 0.1;
    var smpl_m3, smpl_m4 = 0.1;
    var smpl_m3, smpl_m5 = 0.1;
    var smpl_m4, smpl_m5 = 0.1;
end;


// Do not check steady state (un-comment for speed once you know steady state is computed correctly)
options_.steadystate.nocheck = 1;

// Check steady state (un-comment to check whether steady state from .m file satisfies equations.mod)
//steady;
//resid(1)


// Check dynamic regularity conditions (comment out for speed)
//check;
//model_diagnostics;
//model_info;


// Simulate
stoch_simul(irf=0, order=1);

