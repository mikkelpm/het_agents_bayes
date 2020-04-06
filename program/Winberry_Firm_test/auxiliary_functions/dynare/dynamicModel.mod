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
    var smpl_m1_shock; stderr sqrt(0.0548/100); // NEED to change according to cross-sec size, steady state
    var smpl_m2_shock; stderr sqrt(0.0526/100);
    var smpl_m3_shock; stderr sqrt(2*99)*0.0548/100;
    var smpl_m4_shock; stderr sqrt(99*(0.0472^2+0.0548*0.0526))/100;
    var smpl_m5_shock; stderr sqrt(2*99)*0.0526/100;
    var smpl_m1_shock,smpl_m2_shock = 0.0472/100;
    var smpl_m3_shock,smpl_m4_shock = 0.005169552591657*99/100^2; // https://pdfs.semanticscholar.org/ac51/ee74af59c432d493da98bd950cc6f856a0ca.pdf(p4)
    var smpl_m3_shock,smpl_m5_shock = 0.004455426459309*99/100^2;
    var smpl_m4_shock,smpl_m5_shock = 0.004965820793123*99/100^2;
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

// Estimation
estimated_params;
rrhoProd, beta_pdf, 0.5, 1/4;
ssigmaProd, inv_gamma_pdf, 0.04, 1;
end;

estimation(order=1,datafile=simul_moments_N100,mh_replic=10000,mh_jscale=.8);
