function [ys,check] = example1_new_steadystate(ys,exo)
% function [ys,check] = NK_baseline_steadystate(ys,exo)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values fpr the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impos restriction on parameters)

global M_ oo_ 
% warning(['ss ' num2str(oo_.steady_state')])
% disp(['ss ' num2str(oo_.steady_state')])

alpha = 0; % Need to declare alpha here, since Matlab mistakenly thinks alpha is a function otherwise

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

% read out previous steady state (use as initial values)
NumberOfEndogenousVariables = M_.orig_endo_nbr;
for ii = 1:NumberOfEndogenousVariables
  varname = strcat(deblank(M_.endo_names(ii,:)), '_steady');
  eval([ varname ' = oo_.steady_state(' int2str(ii) ');']);
end

% initialize indicator
check = 1;


%% Enter model equations here

optim_options = optimoptions(@fsolve, 'Display', 'off'); % set options for numerical solver

% y and c as fnuctions of k and h
y_fct = @(k,h) (k^alpha)*(h^(1-alpha));
c_fct = @(k,h) (1-alpha)*y_fct(k,h)/(theta*h^(1+psi));

% Function that must be zero at steady state
solve_fct = @(x) [x(1)-beta*(alpha*y_fct(x(1),x(2))+(1-delta)*x(1)); delta*x(1)-(y_fct(x(1),x(2))-c_fct(x(1),x(2)))]; % x(1) is k, x(2) is h

% Solve numerically for k and h
[solve_x,~,exitflag] = fsolve(solve_fct,[k_steady,h_steady],optim_options);
check = (exitflag<1)+0; % =0 if solver converged

% Set new steady state
k = solve_x(1);
h = solve_x(2);
c = c_fct(k,h);
y = y_fct(k,h);
a = 0;
b = 0;

%% end own model equations

% Supply new steady state to Dynare
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end

