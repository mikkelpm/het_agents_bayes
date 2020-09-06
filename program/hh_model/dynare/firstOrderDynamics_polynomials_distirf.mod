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
// MPM: Observables
//----------------------------------------------------------------

varobs logAggregateOutput
@#for iEpsilon in 1 : nEpsilon
    @#for iMoment in 1 : nMeasure
        smpl_m@{iEpsilon}@{iMoment}
    @#endfor
@#endfor
; // order matters

//----------------------------------------------------------------
// 4. Computation
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
steady(nocheck); // Modified wrt firstOrderDynamics_polynomials.mod

// Check regularity conditions (turn on to check)
//check;
//model_diagnostics;
//model_info;

// Simulate 
// Modified wrt firstOrderDynamics_polynomials.mod
stoch_simul(irf=8, order=1, irf_shocks = (aggregateTFPShock), nograph, noprint) 
	lag_moment_1_1 lag_moment_1_2 lag_moment_1_3 lag_moment_2_1 lag_moment_2_2 lag_moment_2_3
	measureCoefficient_1_1 measureCoefficient_1_2 measureCoefficient_1_3 
	measureCoefficient_2_1 measureCoefficient_2_2 measureCoefficient_2_3;
