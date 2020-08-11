% Compute and save steady state

tStart = tic;
fprintf('\nComputing steady state...\n');

check = 0;

%----------------------------------------------------------------
% Read in parameters from Dynare declaration
%----------------------------------------------------------------

% Read out parameters to access them with their name
for iParameter = 1:M_.param_nbr
  paramname = deblank(M_.param_names{iParameter});
%   eval(['global ' paramname]);
  eval([ paramname ' = M_.params(' int2str(iParameter) ');']);
end

%----------------------------------------------------------------
% Steady State
%----------------------------------------------------------------

coreSteadyState;
if check
    return;
end

% Prices
r = aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

%----------------------------------------------------------------
% Save values of steady state variables for Dynare (must be exactly
% as declared in Dynare)
%----------------------------------------------------------------

% Coefficients on conditional expectation function
for iEpsilon = 1 : nEpsilon
    for iAsset = 1 : nAssets
        eval(sprintf('expectationCoefficient_%d_%d = mCoefficients(iEpsilon,iAsset);',...
            iEpsilon,iAsset));
    end
end

% Sample moments of income
mom_lambda = exp((1:nMeasure)*mu_l + 0.5*(1:nMeasure).^2*(-2*mu_l)); % Moments E[lambda^j]
smpl_m = zeros(nEpsilon,nMeasure);
smpl_m(:,1) = w * (mmu * (1-vEpsilonGrid) + (1-ttau)*vEpsilonGrid) + r * mMoments(:,1); % First moment
smpl_m(:,2) = r^2 * mMoments(:,2) * mom_lambda(2) + smpl_m(:,1).^2 * (mom_lambda(2)-1); % Second central moment
smpl_m(:,3) = r^3 * mMoments(:,3) * mom_lambda(3) + 3 * r^2 * mMoments(:,2) .* smpl_m(:,1) * (mom_lambda(3) - mom_lambda(2)) + smpl_m(:,1).^3 * (mom_lambda(3) - 3*mom_lambda(2) + 2); % Third central moment

% Moments, lagged moments, and parameters of density away from borrowing constraint
for iEpsilon = 1 : nEpsilon
    for iMoment = 1 : nMeasure
        eval(sprintf('moment_%d_%d = mMoments(iEpsilon,iMoment);',iEpsilon,iMoment));
        eval(sprintf('lag_moment_%d_%d = moment_%d_%d;',iEpsilon,iMoment,iEpsilon,iMoment));
        eval(sprintf('smpl_m%d%d = smpl_m(iEpsilon,iMoment);',iEpsilon,iMoment,iEpsilon,iMoment));
        eval(sprintf('measureCoefficient_%d_%d = mParameters(iEpsilon,iMoment+1);',iEpsilon,iMoment));
    end
end

% Mass at borrowing constraint
for iEpsilon = 1 : nEpsilon
    eval(sprintf('mHat_%d = mHat(iEpsilon);',iEpsilon));
    eval(sprintf('lag_mHat_%d = mHat_%d;',iEpsilon,iEpsilon));
end

% Other variables
aggregateCapital = (1 - mHat(1,1)) * (1 - aggEmployment) * mMoments(1,1) + (1 - mHat(2,1)) * aggEmployment * mMoments(2,1);
aggregateTFP = 0;
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));
logAggregateInvestment = log(ddelta * aggregateCapital);
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));
logWage = log(w);

% Save steady state values to M_ struct
M_.steady_vars = struct;
for iVar = 1:length(M_.endo_names)
    M_.steady_vars.(M_.endo_names{iVar}) = eval(M_.endo_names{iVar});
end

fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart));
