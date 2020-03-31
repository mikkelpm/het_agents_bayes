% Compute and save steady state 

global M_;

tStart = tic;
fprintf('\nComputing steady state...\n')

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
% Solve for steady state wage
%----------------------------------------------------------------

coreSteadyState;
if check
    return;
end


%----------------------------------------------------------------
% Save grids (may have changed if changed the parameter values)
%----------------------------------------------------------------

% Value function grid
for iState = 1 : nState
    eval(sprintf('valueGrid_%d_1 = mStateGrid(iState,1);',iState));
    eval(sprintf('valueGrid_%d_2 = mStateGrid(iState,2);',iState));
end

% Idiosyncratic shocks
for iShock = 1 : nShocks
	eval(sprintf('shocksGrid_%d = vShocksGrid(iShock);',iShock));
	eval(sprintf('shocksWeights_%d = vShocksWeights(iShock);',iShock));
end

% Quadrature grid
for iState = 1 : nStateQuadrature
	eval(sprintf('quadratureGrid_%d_1 = mQuadratureGrid(iState,1);',iState));
	eval(sprintf('quadratureGrid_%d_2 = mQuadratureGrid(iState,2);',iState));
	eval(sprintf('quadratureWeights_%d = vQuadratureWeights(iState);',iState));
	for iShock = 1 : nShocks
		eval(sprintf('quadratureProdPrimeGrid_%d_%d = mProdPrimeQuadrature(iShock,iState);',iState,iShock));
	end
end

% Value function polynomials
for iState = 1 : nState
	for iPower = 1 : nState
		eval(sprintf('valueFunctionPolys_%d_%d = mStatePoly(iState,iPower);',iState,iPower));
	end
	for iShock = 1 : nShocks
		for iPower = 1 : nProd
			eval(sprintf('valueFunctionPrimePolys_%d_%d_%d = aProdPrimePoly(iShock,iState,iPower);',iState,iShock,iPower));
		end
	end
	for iPower = 1 : nState - nProd
		eval(sprintf('marginalValueFunctionPolys_%d_%d = mStatePoly(iState,iPower);',iState,iPower));
	end
	eval(sprintf('valueFunctionPolySquared_%d = vStatePolySquared(iState);',iState));
end

% Quadrature polynomials
for iState = 1 : nStateQuadrature
	for iPower = 1 : nState
		eval(sprintf('quadraturePolys_%d_%d = mQuadraturePoly(iState,iPower);',iState,iPower));
	end
	for iShock = 1 : nShocks
		for iPower = 1 : nProd
			eval(sprintf('quadraturePrimePolys_%d_%d_%d = aProdPrimeQuadraturePoly(iShock,iState,iPower);',iState,iShock,iPower));
		end
	end
end	


%----------------------------------------------------------------
% Save results for Dynare
%----------------------------------------------------------------

% Compute necessary objects
[vCoefficientsNew,vCapitalAdjust] = updateCoefficients(vCoefficients);
[resid,vMoments,vParameters,aggregateConsumption,marginalUtility,aggregateOutput,aggregateCapital,aggregateInvestment] = ...
	computeLMCResidualPolynomials(wage,vParameters,vMomentsHistogram,mGridMoments);
aggregateHours = resid + nSS;
	
% Value function
for iState = 1 : nState
	eval(sprintf('valueCoefficient_%d = marginalUtility * vCoefficientsNew(iState);',iState));
end

% Adjust capital policy function
for iState = 1 : nState
	eval(sprintf('capitalAdjust_%d = vCapitalAdjust(iState);',iState));
end

% Moments
for iMoment = 1 : nMeasureCoefficients
	eval(sprintf('moment_%d = vMoments(iMoment);',iMoment));
end

% NEED to change for different nMeasure
smpl_m1 = -(moment_1-log(wage)+log(nnu)+ttheta*moment_2)/(nnu-1);
smpl_m2 = moment_2;
smpl_m3 = (moment_3+ttheta^2*moment_5+2*ttheta*moment_4)/(nnu-1)^2;
smpl_m4 = -(ttheta*moment_5+moment_4)/(nnu-1);
smpl_m5 = moment_5;

% Parameters
for iParameter = 0 : nMeasureCoefficients
	eval(sprintf('measureCoefficient_%d = vParameters(iParameter+1);',iParameter));
end

% Aggregate shocks
aggregateTFP    	= 0;
aggregateQ      	= 0;
smpl_m1_shock = 0;
smpl_m2_shock = 0;
smpl_m3_shock = 0;
smpl_m4_shock = 0;
smpl_m5_shock = 0;

% Other auxiliary variables of interest
expectedMarginalUtilityPrime 	= marginalUtility;
realInterestRate 					= 1 / bbeta;
logAggregateOutput 				= log(aggregateOutput);
logAggregateInvestment 		= log(aggregateInvestment);
logAggregateHours 				= log(aggregateHours);
logAggregateConsumption 		= log(aggregateConsumption);
logWage 								= log(wage);

vCoefficients 							= marginalUtility * vCoefficientsNew;
vParameters 							= vParameters(2:nMeasureCoefficients+1);

wageSS 								= wage;
marginalUtilitySS 					= marginalUtility;
cchi 										= (wage * marginalUtility / (aggregateHours ^ pphi));

logMarginalUtility					= log(marginalUtilitySS);

logAggregateOutputObs = logAggregateOutput;
logAggregateInvestmentObs = logAggregateInvestment;

% Save
save_vars = cellstr(M_.endo_names);
save('steady_vars.mat', save_vars{:}, 'save_vars');

% Update parameters which were changed in the steady state file
for iParameter = 18:length(M_.params)                               % start at 18 ( = cchi) so don't reset other parameters, in case those are reset elsewhere
  eval([ 'M_.params(' num2str(iParameter) ') = ' M_.param_names(iParameter,:) ';' ])
end

% Print elapsed time
fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart))