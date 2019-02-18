function [ys,check] = firstOrderDynamics_polynomials_steadystate(ys,exo)

% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
%
% Thomas Winberry, July 26th, 2016

tStart = tic;
fprintf('\nComputing steady state...\n')

%----------------------------------------------------------------
% Call parameters (the next set of commands will overwrite some)
%----------------------------------------------------------------
load('economicParameters.mat')
load('approximationParameters.mat')

%----------------------------------------------------------------
% Read in parameters from Dynare declaration
%----------------------------------------------------------------

% Initialize indicator
check = 0;

% Read parameters from Dynare
global M_ % LL: seems that we need to keep this global to satisfies the requirement of Dynare...

% Read out parameters to access them with their name
for iParameter = 1:M_.param_nbr-3
  paramname = deblank(M_.param_names(iParameter,:));
%   eval(['global ' paramname]);
  eval([ paramname ' = M_.params(' int2str(iParameter) ');']);
end

%----------------------------------------------------------------
% Steady State
%----------------------------------------------------------------
displayOpt = 'off';       % 'iter-detailed' or 'off'
coreSteadyState;

w = w_SS;
A = A_SS;
r = r_SS;
i = r_SS;
N = N_SS;
d = d_SS;
T = T_SS;
ppi = 0;

M_.params(M_.param_nbr-2) = r_SS;
M_.params(M_.param_nbr-1) = G_SS;
M_.params(M_.param_nbr) = B_SS;

%----------------------------------------------------------------
% Save values of steady state variables for Dynare (must be exactly
% as declared in Dynare)
%----------------------------------------------------------------

% Coefficients on conditional expectation function
for iShare = 1:nShare
    for iz = 1 : nz
        for iAsset = 1 : nAssets
            eval(sprintf('expectationCoefficient_%d_%d_%d = mCoefficients(iz,iAsset,iShare);',...
                iz,iAsset,iShare));
        end
    end
end

% Moments and parameters of density away from borrowing constraint
for iShare = 1:nShare
    for iz = 1 : nz
        for iMoment = 1 : nMeasure
            eval(sprintf('moment_%d_%d_%d = mMoments(iz,iMoment,iShare);',iz,iMoment,iShare));
            eval(sprintf('measureCoefficient_%d_%d_%d = mParameters(iz,iMoment+1,iShare);',iz,iMoment,iShare));
        end
    end
end

% Mass at borrowing constraint
for iShare = 1:nShare
    for iz = 1 : nz
        eval(sprintf('mHat_%d_%d = mHat(iz,iShare);',iz,iShare));
    end
end

% Save endogenous variables back into ys
for ii = 1 : M_.orig_endo_nbr;
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']); 
end

fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart))
