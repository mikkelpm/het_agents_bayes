function [residual,vCoefficientsOptional,vCapitalCoefficientsOptional,vCapitalAdjustOptional,...
	vCapitalConstrainedOptional,vCutoffOptional,vHistogramOptional] = computeLMCResidualHistogram(wage);

% Computes the residual of the labor market clearing condition, approximating the distribution
% using a histogram; used to compute initial guess of parametric family
%
% Inputs
%	(1) wage: candidate wage
%
% Outputs
%	(1) residual: residual of market clearing condition (to use in root finder)
%	(2) (optional) CoefficientsOptional: coefficients on value function 
%	(3) (optional) vCapitalCoefficientsOptional: coefficients on adjust capital decision 
%	(4) (optional) vCapitalAdjustOptional: capital decision if adjust, along histogram grid 
%	(5) (optional) vCapitalConstrainedOptional: capital decision if not adjust, along histogram grid 
%	(6) (optional) vCutoffOptional: adjustment cutoff, along histogram grid
%
% Thomas Winberry, February 14th, 2018

% Declare global variables used in this function
global mStateGrid ttheta nnu ddelta mStatePoly vStatePolySquared tolerance maxIterations nState mFinePoly mFineGrid ...
	aaUpper aaLower capitalMin capitalMax nShocks aProdPrimeFinePoly nStateFine nProd nCapital bbeta dampening ...
	nProdFine nCapitalFine vShocksWeights ppsiCapital nSS
	
%----------------------------------------------------------------
% Compute value function
%----------------------------------------------------------------

% Compute labor demand over the grid (independent of value function, but depends on wage)
global vLaborDemandGrid vProfitGrid
vLaborDemandGrid 	= ((exp(mStateGrid(:,1)) .* (mStateGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu));
vProfitGrid 				= exp(mStateGrid(:,1)) .* (mStateGrid(:,2) .^ ttheta) .* (vLaborDemandGrid .^ nnu) - wage * vLaborDemandGrid;

% Initialize value function
vInitGrid 				= vProfitGrid + (1 - ddelta) * mStateGrid(:,2);
vCoefficients 			= sum(mStatePoly' .* (ones(nState,1) * vInitGrid'),2);
vCoefficients 			= vCoefficients ./ vStatePolySquared;

% Initialize iteration
err 						= 100; 
iterations 				= 1;
t0 						= tic;

% Do value function iteration
while err > tolerance && iterations <= maxIterations

	[vCoefficientsNew,vCapitalAdjust] = updateCoefficients(vCoefficients);
	err 											= max(abs(vCoefficientsNew - vCoefficients));
	iterations 									= iterations + 1;
	vCoefficients 								= dampening * vCoefficients + (1 - dampening) * vCoefficientsNew;

end


%----------------------------------------------------------------
% Compute value and policy functions over fine grid
%----------------------------------------------------------------

% Compute polynomial approximation of capital accumulation policy conditional on adjustment
vCapitalCoefficients = sum(mStatePoly' .* (ones(nState,1) * vCapitalAdjust'),2);
vCapitalCoefficients = vCapitalCoefficients ./ vStatePolySquared;

% Compute policy functions over fine grid
[vCapitalAdjust,vCapitalConstrained,vCutoff] = computePolicies(vCoefficients,vCapitalCoefficients,wage,mFineGrid,mFinePoly,aProdPrimeFinePoly);


%----------------------------------------------------------------
% Compute stationary distribution from these decision rules,
% using histogram 
%----------------------------------------------------------------

% Compute discrete transition matrix
mTransition 	= sparse(computeDiscreteTransitionMatrix(vCapitalAdjust,vCapitalConstrained,vCutoff));

% Iterate over transition matrix
err 				= 100;
iterations 		= 0;
vHistogram 	= ones(nStateFine,1) ./ nStateFine;
while err > 1e-16 && iterations < 1e4

    vHistogramNew 	= mTransition' * vHistogram;
    err 					= max(abs(vHistogramNew - vHistogram));
    iterations 			= iterations + 1;
    vHistogram 		= vHistogramNew;
    
end

%----------------------------------------------------------------
% Compute output of the function
%----------------------------------------------------------------

% Market clearing residual
vLaborDemand 			= ((exp(mFineGrid(:,1)) .* (mFineGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu)) + ((vCutoff .^ 2) ./ (2 * ppsiCapital));
residual 					= sum(vLaborDemand .* vHistogram) - nSS;		% will later calibrate cchi to ensure labor supply = nSS

% Optional output
if nargout > 1
	
	vCoefficientsOptional 			= vCoefficients;
	vCapitalCoefficientsOptional = vCapitalCoefficients;
	vCapitalAdjustOptional 		= vCapitalAdjust;
	vCapitalConstrainedOptional = vCapitalConstrained;
	vCutoffOptional 					= vCutoff;
	vHistogramOptional 			= vHistogram;

end	
