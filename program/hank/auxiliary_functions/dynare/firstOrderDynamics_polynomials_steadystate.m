function [ys,check] = firstOrderDynamics_polynomials_steadystate(ys,exo)

% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
%
% Thomas Winberry, July 26th, 2016

% Load steady state
try
    load('steady_vars.mat');
    check = 0;
catch
    computeSteadyState;
end

% Save endogenous variables back into ys
for ii = 1 : length(save_vars) % M_.orig_endo_nbr
  varname = deblank(save_vars{ii}); % %M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']); 
end
