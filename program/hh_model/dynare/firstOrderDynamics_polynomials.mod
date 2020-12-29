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
// Observables
//----------------------------------------------------------------

varobs logAggregateOutput
@#for iEpsilon in 1 : nEpsilon
    @#for iMoment in 1 : nMeasure
        smpl_m@{iEpsilon}@{iMoment}
    @#endfor
@#endfor
; // order matters

//----------------------------------------------------------------
// Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
    var aggregateTFPShock = 1;

    // Placeholder values for measurement error
    var logAggregateOutput = 1;
    @#for iMoment in 1 : nMeasure
        var smpl_m1@{iMoment} = 1;
        var smpl_m2@{iMoment} = 1;
    @#endfor
	@#for iMoment in 1 : nMeasure
        @#for jMoment in iMoment+1 : nMeasure
			var smpl_m1@{iMoment},smpl_m1@{jMoment} = .1;
			var smpl_m2@{iMoment},smpl_m2@{jMoment} = .1;
		@#endfor
    @#endfor
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
