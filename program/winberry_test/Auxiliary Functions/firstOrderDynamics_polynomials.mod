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
var logAggregateOutputObs logWageObs;
varexo measErr measErr_wage;

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_polynomials.mod"
logAggregateOutputObs = logAggregateOutput + measErr;
logWageObs = logWage + measErr_wage;

end;

//----------------------------------------------------------------
// MPM: Observables
//----------------------------------------------------------------

//varobs logAggregateOutputObs;
varobs logAggregateOutputObs logWageObs;

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
    var aggregateTFPShock = 1;
    var measErr = 0.03^2;
    var measErr_wage = 0.03^2;
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
