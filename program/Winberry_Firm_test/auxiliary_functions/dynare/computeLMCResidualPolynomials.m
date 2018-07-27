function [residual,vMoments,vParameters,aggregateConsumption,marginalUtility,aggregateOutput,aggregateCapital,aggregateInvestment] = ...
	computeLMCResidualPolynomials(wage,vParameters,vMoments,mGridMoments);

% Computes the residual of the labor market clearing condition, approximating the distribution
% using exponential polynomials
% 
% Inputs
%   (1) wage: candidate wage
%	(2) vParmaters: initial guess of distribution parameters
%	(3) vMoments: initial guess of distribution moments
%	(4) mGridMoments: grid of centralized moments, corresponding to vMoments
%
% Outputs
%   (1) resid: aggregate labor demand - supply
%	(2) vMoments: moments of parametric family (optional)
%	(3) vParameters: parameters of parametric family (optional)
%	(4) aggregateConsumption (optional)
%	(5) marginalUtility (optional)
%	(6) aggregateOutput (optional)
%	(7) aggregateCapital (optional)
%	(8) aggregateInvestment (optional)
% 
% Thomas Winberry, September 8th, 2016

% Declare global variables
global mStateGrid ttheta nnu ddelta mStatePoly vStatePolySquared tolerance maxIterations nState mFinePoly mFineGrid ...
	aaUpper aaLower capitalMin capitalMax nShocks aProdPrimeFinePoly nStateFine nProd nCapital bbeta ppsi2Capital dampening ...
	nProdFine nCapitalFine vShocksWeights ppsiCapital nProdQuadrature nCapitalQuadrature nStateQuadrature mQuadratureGrid ...
	mQuadraturePoly aProdPrimeQuadraturePoly nMeasure nMeasureCoefficients mProdPrimeQuadrature vQuadratureWeights nSS ssigma

	
%----------------------------------------------------------------
% Compute value function
%----------------------------------------------------------------

global vLaborDemandGrid vProfitGrid
vLaborDemandGrid = ((exp(mStateGrid(:,1)) .* (mStateGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu));
vProfitGrid = exp(mStateGrid(:,1)) .* (mStateGrid(:,2) .^ ttheta) .* (vLaborDemandGrid .^ nnu) - wage * vLaborDemandGrid;

% Initialize
vInitGrid = vProfitGrid + (1 - ddelta) * mStateGrid(:,2);
vCoefficients = sum(mStatePoly' .* (ones(nState,1) * vInitGrid'),2);
vCoefficients = vCoefficients ./ vStatePolySquared;
err = 100; iterations = 1;
t0 = tic;

% Iterate
while err > tolerance && iterations <= maxIterations

	[vCoefficientsNew,vCapitalAdjust] = updateCoefficients(vCoefficients);
	err = max(abs(vCoefficientsNew - vCoefficients));
	iterations = iterations + 1;
	vCoefficients = dampening * vCoefficients + (1 - dampening) * vCoefficientsNew;

end


%----------------------------------------------------------------
% Compute value and policy functions over quadrature grid
%----------------------------------------------------------------

% Compute polynomial approximation of investment policy conditional on adjustment
vCapitalCoefficients = sum(mStatePoly' .* (ones(nState,1) * vCapitalAdjust'),2);
vCapitalCoefficients = vCapitalCoefficients ./ vStatePolySquared;

% Compute policy functions
[vCapitalAdjust,vCapitalConstrained,vCutoff] = computePolicies(vCoefficients,vCapitalCoefficients,wage,mQuadratureGrid,mQuadraturePoly,aProdPrimeQuadraturePoly);

% Enforce choices to be within bounds
vCapitalAdjust = min(max(capitalMin * ones(nStateQuadrature,1),vCapitalAdjust),capitalMax * ones(nStateQuadrature,1));
vCapitalConstrained = min(max(capitalMin * ones(nStateQuadrature,1),vCapitalConstrained),capitalMax * ones(nStateQuadrature,1));

% Replicate choices for each draw of future idioysncratic shock
vCapitalAdjustPrime = repmat(vCapitalAdjust', [nShocks 1]);
vCapitalConstrainedPrime = repmat(vCapitalConstrained', [nShocks 1]);
vCutoffPrime = repmat(vCutoff', [nShocks 1]);


%----------------------------------------------------------------
% Compute stationary distribution by iterating on law of motion
%----------------------------------------------------------------

% Initialize objects for the iteration
err = 100; iteration = 1; maxIterations = 200; dampening = 0;
options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
	'MaxFunEvals',50000,'TolFun',1e-6,'GradObj','on','MaxIter',1000);

% Iterate
while err > 1e-6 && iteration <= maxIterations

	% Compute parameters of distribution
	objectiveFunction 				= @(vParametersTilde) parametersResidual(vParametersTilde,mGridMoments);
	[vParameters,normalization] = fminunc(objectiveFunction,vParameters(2:nMeasureCoefficients+1,1),options);
    %[vParameters,normalization] = fminunc(objectiveFunction,zeros(nMeasureCoefficients,1),options);		% if above initial guess no good, try this one
	vParameters 						= [1 / normalization; vParameters];

	% Compute new moments and centered moments grid
	% First moments
	vMomentsNew = zeros(2,1);
	for i = 0:1
		mIntegrand 				= (mProdPrimeQuadrature .^ (1 - i)) .* ((vCutoffPrime ./ ppsiCapital) .* (log(vCapitalAdjustPrime) .^ i) + ...
											(1 - (vCutoffPrime ./ ppsiCapital)) .* (log(vCapitalConstrainedPrime) .^ i));
		vIntegrand 				= (vShocksWeights' * mIntegrand)';
		vMomentsNew(i+1,1) 	= vQuadratureWeights' * (vIntegrand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)));
	end
	mGridMomentsNew 			= [mQuadratureGrid(:,1) - vMomentsNew(1,1) log(mQuadratureGrid(:,2)) - vMomentsNew(2,1)];

	% Higher-order moments
	for i = 2:nMeasure
		for j = 0:i
			mIntegrand 			= (mProdPrimeQuadrature .^ (i - j)) .* ((vCutoffPrime ./ ppsiCapital) .* ((log(vCapitalAdjustPrime) - vMomentsNew(2,1)) .^ j) + ...
											(1 - (vCutoffPrime ./ ppsiCapital)) .* ((log(vCapitalConstrainedPrime) - vMomentsNew(2,1)) .^ j));
			vIntegrand 			= (vShocksWeights' * mIntegrand)';
			vMomentsNew 		= [vMomentsNew; vQuadratureWeights' * (vIntegrand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)))];
			mGridMomentsNew 	= [mGridMomentsNew (((mQuadratureGrid(:,1) - vMomentsNew(1,1)) .^ (i - j)) .* ((log(mQuadratureGrid(:,2)) - ...
											vMomentsNew(2,1)) .^ j) - vQuadratureWeights' * (vIntegrand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1))))];
		end
	end
		
	% Update iteration
	err 					= max(abs(vMomentsNew - vMoments)); 
	iteration 			= iteration + 1;
	vMoments 			= (1 - dampening) * vMomentsNew + dampening * vMoments;
	mGridMoments 	= (1 - dampening) * mGridMomentsNew + dampening * mGridMoments;
	
end


%----------------------------------------------------------------
% Compute labor market clearing residual
%----------------------------------------------------------------

% Compute labor demand along grid
vLaborDemand 		= ((exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu)) + ...
								((vCutoff .^ 2) ./ (2 * ppsiCapital));

% Integrate
aggLaborDemand 	= vQuadratureWeights' * (vLaborDemand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)));
residual 				= aggLaborDemand - nSS;


%----------------------------------------------------------------
% Optional: compute steady state aggregates
%----------------------------------------------------------------

if nargout > 1

	% Aggregate consumption 
	vLaborDemand 			= ((exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) * nnu) / wage) .^ (1 / (1 - nnu));
	vIntegrand 				= exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) .* (vLaborDemand .^ nnu) + ...
										(vCutoff ./ ppsiCapital) .* (-(vCapitalAdjust - (1 - ddelta) * mQuadratureGrid(:,2))) + (1 - (vCutoff ./ ppsiCapital)) .* ...
										(-(vCapitalConstrained - (1 - ddelta) * mQuadratureGrid(:,2)));
	aggregateConsumption = vQuadratureWeights' * (vIntegrand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)));

	% Marginal utility
	marginalUtility 			= aggregateConsumption ^ (-ssigma);

	% Output
	vIntegrand 				= exp(mQuadratureGrid(:,1)) .* (mQuadratureGrid(:,2) .^ ttheta) .* (vLaborDemand .^ nnu);
	aggregateOutput 		= vQuadratureWeights' * (vIntegrand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)));

	% Capital
	aggregateCapital 		= vQuadratureWeights' * (mQuadratureGrid(:,2) .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)));

	% Investment
	vIntegrand 				= (vCutoff ./ ppsiCapital) .* (vCapitalAdjust - (1 - ddelta) * mQuadratureGrid(:,2)) + (1 - (vCutoff ./ ppsiCapital)) .* ...
										(vCapitalConstrained - (1 - ddelta) * mQuadratureGrid(:,2));
	aggregateInvestment 	= vQuadratureWeights' * (vIntegrand .* vParameters(1,1) .* exp(mGridMoments * vParameters(2:nMeasureCoefficients+1)));	

end
