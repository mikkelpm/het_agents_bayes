% Computes market clearing wage and associated distribution and decision rules
%
% Thomas Winberry, February 14th, 2018

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
load_mat('grids');

% Polynomials 
load_mat('polynomials');


%----------------------------------------------------------------
% Compute initial guess of wage from histogram
% (will use as initial guess for parametric distribution below)
%----------------------------------------------------------------

global wage
wage 	= wRepSS;			% Initial guess from representative agent model

t0 		= tic;
fprintf('Computing initial guess from histogram...\n')

% Solve for market clearing wage
f 								= @(wage) computeLMCResidualHistogram(wage);
options 						= optimoptions('fsolve','Display','off');		% set 'off' if don't want to print progress
[wageInit,err,exitflag,output] 	= fsolve(f,wage,options);

% Compute decisions and distribution at wageInit
[~,vCoefficients,vCapitalCoefficients,vCapitalAdjust,...
	vCapitalConstrained,vCutoff,vHistogram] 	= computeLMCResidualHistogram(wageInit);

	
%----------------------------------------------------------------
% Compute moments of histogram to use as initial guess for parametric family
%----------------------------------------------------------------

% Compute moments from histogram and use to create grid around centered moments
vMomentsHistogram 	    = zeros(2,1);
vMomentsHistogram(1,1) = sum(mFineGrid(:,1) .* vHistogram);			% first moment of productivity
vMomentsHistogram(2,1) = sum(log(mFineGrid(:,2)) .* vHistogram);		% first moment of log capital

mGridMoments 				= zeros(nStateQuadrature,2);
mGridMoments(:,1) 			= mQuadratureGrid(:,1) - vMomentsHistogram(1,1);
mGridMoments(:,2) 			= log(mQuadratureGrid(:,2)) - vMomentsHistogram(2,1);

for i = 2:nMeasure		%	higher-order moments
	for j = 0:i	
		newMoment 			= sum(((mFineGrid(:,1) - vMomentsHistogram(1,1)) .^ (i - j)) .* ...
											((log(mFineGrid(:,2)) - vMomentsHistogram(2,1)) .^ j) .* vHistogram);
		vMomentsHistogram= [vMomentsHistogram; newMoment];
		mGridMoments 		= [mGridMoments ((mQuadratureGrid(:,1) - vMomentsHistogram(1,1)) .^ (i - j)) .* ...
											((log(mQuadratureGrid(:,2)) - vMomentsHistogram(2,1)) .^ j) - newMoment];
	end
end

% Compute parameters by minimizing integral equation
options 								= optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
											'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
objectiveFunction 				= @(vParametersTilde) parametersResidual(vParametersTilde,mGridMoments);
[vParameters,normalization] = fminunc(objectiveFunction,zeros(nMeasureCoefficients,1),options);
%[vParameters,normalization] = fminunc(objectiveFunction,-1e-2*ones(nMeasureCoefficients,1),options);  % if optimizer is not converging or has problems with that initial guess, try this initial guess
vParameters 						= [1 / normalization; vParameters];

fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))


%----------------------------------------------------------------
% Compute refined wage using exponential polynomials
%----------------------------------------------------------------

t0 = tic;
fprintf('Compute refined steady state from exponential polynomials...\n')

% Solve for market clearing wage
f 					= @(wage) computeLMCResidualPolynomials(wage,vParameters,vMomentsHistogram,mGridMoments);
options 			= optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
if abs(f(wageInit)) > 1e-4
	[wage,err,exitflag] = fsolve(f,wageInit,options);
end

% Return error if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end	

fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
