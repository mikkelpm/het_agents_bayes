// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating conditional expectation with 
// polynomials)
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters_polynomials.mod"

//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables_polynomials.mod"
// MPM: include observation with measurement error
var logAggregateOutputObs;
varexo measErr;

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_polynomials.mod"
logAggregateOutputObs = logAggregateOutput + ssigmaMeas*measErr; // MPM: measurement error

end;

//----------------------------------------------------------------
// MPM: Observables
//----------------------------------------------------------------

varobs logAggregateOutputObs
@#for iMoment in 1 : nMeasure
    smpl_m1@{iMoment} smpl_m2@{iMoment}
@#endfor
;

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
    var aggregateTFPShock = 1;
    // MPM: Measurement error
    var measErr = 1;

    // Placeholder values for measurement error
    @#for iMoment in 1 : nMeasure
        var smpl_m1@{iMoment} = 1;
        var smpl_m2@{iMoment} = 1;
    @#endfor
    // OK for nMeasure = 2 or 3, NEED to change if nMeasure takes other values
    @#if nMeasure <= 3
        var smpl_m11, smpl_m13 = 0.1;
        var smpl_m21, smpl_m23 = 0.1;
    @#endif
end;

options_.steadystate.nocheck = 1;

// Compute steady state (nocheck option ensures that Dynare runs even if steady
// state only computed approximately, i.e., with small numerical error)
//steady(nocheck);

// Check regularity conditions (turn on to check)
//check;
//model_diagnostics;
//model_info;

// Simulate
stoch_simul(irf=0, order=1);
