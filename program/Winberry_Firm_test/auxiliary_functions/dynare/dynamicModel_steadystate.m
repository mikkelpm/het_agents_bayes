function [ys,params,check] = dynamicModel_steadystate(ys,exo,M_,options_)

% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
%
% Thomas Winberry, February 15th, 2018

% Load steady state
try
    load('steady_vars.mat');
catch
    setParameters;
    compute_steady_state;
end

params = M_.params;

% Save endogenous variables back into ys
for ii = 1 : length(save_vars) % M_.orig_endo_nbr
  varname = deblank(save_vars{ii}); % %M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']); 
end

check = 0;