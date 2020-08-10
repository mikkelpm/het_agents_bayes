function [ys,params,check] = firstOrderDynamics_polynomials_steadystate(ys,exo,M_,options_)

% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
%
% Thomas Winberry, July 26th, 2016

% Compute steady state if not already computed
if ~isfield(M_, 'steady_vars')
    setParameters;
    compute_steady_state;
end

params = M_.params;

% Save endogenous variables back into ys
for iVar = 1:length(M_.endo_names)
    ys(iVar) = M_.steady_vars.(M_.endo_names{iVar});
end

check = 0;