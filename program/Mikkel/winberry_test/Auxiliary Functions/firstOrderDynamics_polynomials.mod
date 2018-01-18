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
var logAggregateOutputObs;
varexo measErr;

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_polynomials.mod"
logAggregateOutputObs = logAggregateOutput + measErr;

end;

//----------------------------------------------------------------
// MPM: Observables
//----------------------------------------------------------------

varobs logAggregateOutputObs;
//logAggregateOutput logAggregateInvestment logAggregateConsumption logWage

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
    var aggregateTFPShock = 1;
    var measErr = 1;
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
