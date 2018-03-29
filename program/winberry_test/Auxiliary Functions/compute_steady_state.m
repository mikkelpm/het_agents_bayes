% Compute and save steady state

global M_;

tStart = tic;
fprintf('\nComputing steady state...\n');

check = 0;

%----------------------------------------------------------------
% Read in parameters from Dynare declaration
%----------------------------------------------------------------

% Read out parameters to access them with their name
for iParameter = 1:M_.param_nbr
  paramname = deblank(M_.param_names(iParameter,:));
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

% Moments and parameters of density away from borrowing constraint
for iEpsilon = 1 : nEpsilon
    for iMoment = 1 : nMeasure
        eval(sprintf('moment_%d_%d = mMoments(iEpsilon,iMoment);',iEpsilon,iMoment));
        eval(sprintf('measureCoefficient_%d_%d = mParameters(iEpsilon,iMoment+1);',iEpsilon,iMoment));
    end
end

% Mass at borrowing constraint
for iEpsilon = 1 : nEpsilon
    eval(sprintf('mHat_%d = mHat(iEpsilon);',iEpsilon));
end

% Other variables
aggregateCapital = (1 - mHat(1,1)) * (1 - aggEmployment) * mMoments(1,1) + (1 - mHat(2,1)) * aggEmployment * mMoments(2,1);
aggregateTFP = 0;
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));
logAggregateInvestment = log(ddelta * aggregateCapital);
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));
logWage = log(w);
logAggregateOutputObs = logAggregateOutput;

% Save
save_vars = [cellstr(M_.endo_names); {'M_'}];
save('steady_vars.mat', save_vars{:});

fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart));
