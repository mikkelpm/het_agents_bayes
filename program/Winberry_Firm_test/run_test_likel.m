clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');


%% Settings

% Decide what to do
is_data_gen = 0; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation
is_profile = 0; %whether run profiler for execution time
if is_profile == 1
    profile on
end

% Model/data settings
T = 50;                                % Number of periods of simulated macro data
ts_micro = 10:10:T;                        % Time periods where we observe micro data
N_micro = 1e2;                             % Number of micro entities per non-missing time period

% Parameter values to check (prod dynamics)
param1_vals = 0.011*[0.5 1 2]; %.53*[0.9 1 1.1]; %.53*[0.5 1 2];
param2_vals = 0.0083;%*[0.5 1 2]; %.0364*[0.9 1 1.1]; %.0364*[0.5 1 2];

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 201810011;                    % Random number generator seed for initial simulation


%% Set economic parameters 

global ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ  cchi Nmicro

% Technology
ttheta 			= .256;								% capital coefficient
nnu 				= .64;								% labor coefficient
ddelta 			= .085;								% depreciation (annual)
rrhoProd 		= .53; %.859; 								% persistence of idiosyncratic shocks (annual)
ssigmaProd 	= .0364;  %.022;								% SD innovations of idiosycnratic shocks (annual)
aaUpper 		= .011; 								% no fixed cost region upper bound
aaLower 		= -.011;								% no fixed cost region lower bound
ppsiCapital 	= .0083;							% upper bound on fixed adjustment cost draws
% ppsiCapital     = 1e-5; 

% Preferences
bbeta 			= .961;								% discount factor (annual)
ssigma 			= 1;									% coefficient of relative risk aversion
pphi 				= 1 / 1e5;							% inverse Frisch elasticity of labor supply
nSS 				= 1 / 3;								% hours worked in steady state
cchi				= 1;									% labor disutility (will be calibrated in Dynare's steady state file to ensure labor supply = nSS)

% Aggregate shocks
rrhoTFP			= 0.859;							% persistence of aggregate TFP (annual)
ssigmaTFP		= .014;								% SD of innovations of aggregate TFP (annual)
rrhoQ				= 0.859;							% persistence of aggregate investment-specific shock (annual)
ssigmaQ		= .014;								% SD of innovations of aggregate investment-specific shock (annual)
corrTFPQ		= 0;									% loading on TFP shock in evolution of investment-specific shock

%% Set approximation parameters

global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
	maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
	nCapitalQuadrature kRepSS wRepSS

% Order of approximation of value function
nProd 			= 3;										% order of polynomials in productivity
nCapital 		= 5;										% order of polynomials in capital
nState 			= nProd * nCapital;					% total number of coefficients

% Shocks 
nShocks 		= 3;										% order of Gauss-Hermite quadrature over idiosyncratic shocks

% Finer grid for analyzing policy functions and computing histogram
nProdFine 		= 60;
nCapitalFine 	= 40;
nStateFine 	= nProdFine * nCapitalFine;

% Iteration on value function
maxIterations	= 100;
tolerance 		= 1e-6;
acc 				= 500;									% number of iterations in "Howard improvement step"
dampening 	= 0;										% weight on old iteration in updating step

% Approximation of distribution
nMeasure 				= 2;							% order of polynomial approximating distribution
nProdQuadrature 		= 8; 							% number of quadrature points in productivity dimension
nCapitalQuadrature 	= 10;						% number of quadrature points in capital dimension
nStateQuadrature 		= nProdQuadrature * nCapitalQuadrature;
nMeasureCoefficients 	= (nMeasure * (nMeasure + 1)) / 2 + nMeasure;

%% Save parameters

cd('./auxiliary_functions/dynare');
delete steady_vars.mat;
saveParameters;


%% Initial Dynare run

rng(rng_seed);
dynare dynamicModel noclearall nopathchange; % Run Dynare once to process model file
compute_meas_err;
return;


%% Simulate data

if is_data_gen == 0
    
    % Load previous data
    load('simul.mat')
    load('simul_data_micro.mat');
%     load('simul_data_micro_indv_param.mat');
    
else
    
    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw micro data
    simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, num_interp);
    save('simul_data_micro.mat','simul_data_micro');
    
%     % draw individual productivities and incomes
%     simul_data_micro_indv_param = simulate_micro_indv_param(simul_data_micro);
%     save('simul_data_micro_indv_param.mat','simul_data_micro_indv_param');
    
end


%% Compute likelihood

loglikes = nan(length(param1_vals),length(param2_vals));
loglikes_macro = nan(length(param1_vals),length(param2_vals));
loglikes_micro = nan(length(param1_vals),length(param2_vals));

disp('Computing likelihood...');
timer_likelihood = tic;

if contains(pwd,'Laura')
    poolobj = parpool(2);
else
    poolobj = parpool;
end

for iter_i=1:length(param1_vals) % For each parameter...
    
    for iter_j=1:length(param2_vals) % For each parameter...
        
        % Set new parameters
        aaUpper = param1_vals(iter_i);
        aaLower = -aaUpper;
        ppsiCapital = param2_vals(iter_j);
%         rrhoProd = param1_vals(iter_i);
%         ssigmaProd = param2_vals(iter_j);

        fprintf(['%s' repmat('%6.4f ',1,2),'%s\n'], '[aaUpper,ppsiCapital] = [',aaUpper,ppsiCapital,']');
%         fprintf(['%s' repmat('%6.4f ',1,2),'%s\n'], '[rrhoProd,ssigmaProd] = [',rrhoProd,ssigmaProd,']');

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
        compute_Sigma_e;        % Update measurement error var-cov matrix

        % Log likelihood of proposal
        [loglikes(iter_i,iter_j), loglikes_macro(iter_i,iter_j), loglikes_micro(iter_i,iter_j)] = ...
            loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                                       num_smooth_draws, num_interp, num_burnin_periods, ...
                                       M_, oo_, options_);
    
    end
    
end

delete(poolobj);

likelihood_elapsed = toc(timer_likelihood);
fprintf('%s%8.2f\n', 'Done. Elapsed minutes: ', likelihood_elapsed/60);

cd('../../');
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');

if is_profile
    profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);
end
