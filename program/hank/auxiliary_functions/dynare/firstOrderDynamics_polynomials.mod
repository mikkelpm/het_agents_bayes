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

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_polynomials.mod"

end;

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

shocks;
    var epsilonA = 1;
	var epsilonM = 1;
end;

//----------------------------------------------------------------
// Observables
//----------------------------------------------------------------

varobs Y w;

//----------------------------------------------------------------
// Computation
//----------------------------------------------------------------

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
//stoch_simul(irf=40, order=1) A i r w Y N C T d ppi;
