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
var logAggregateOutputObs logAggregateInvestmentObs;
varexo measErrOutput measErrInvestment;


//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

	@#include "equations.mod"
    logAggregateOutputObs = logAggregateOutput + 0.02*measErrOutput;
    logAggregateInvestmentObs = logAggregateInvestment + 0.02*measErrInvestment;

end;

varobs logAggregateOutputObs logAggregateInvestmentObs;

//----------------------------------------------------------------
// Set options
//----------------------------------------------------------------

// Specify shock process (=1 to include shock; =0 to not include shock)
shocks;
    var aggregateTFPShock 	= 1;
    var aggregateQShock   	= 1;
    var measErrOutput = 1;
    var measErrInvestment = 1;
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
