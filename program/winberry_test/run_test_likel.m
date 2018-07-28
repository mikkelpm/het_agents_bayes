clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');


%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation

% Model/data settings
T = 100;                                % Number of periods of simulated macro data
ts_hh = 20:20:T;                        % Time periods where we observe micro data
N_hh = 1e3;                             % Number of households per non-missing time period

% Parameter values to check
param1_vals = [0.93 0.96 0.99];
param2_vals = [-0.5 -0.25 -0.1]; % [0.01 0.02 0.03];

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20180727;                    % Random number generator seed for initial simulation

%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP ttau mu_l ssigmaMeas;
	
% Preferences
bbeta = .96;										% discount factor (annual calibration)
ssigma = 1;											% coefficient of relative risk aversion
aaBar = 0;											% borrowing constraint

% Technology
aalpha = .36;										% capital share
ddelta = .1;										% depreciation rate (annual calibration)

% Idiosyncratic shocks
vEpsilonGrid = [0;1];
aggEmployment = .93;
uDuration = 1;

% Unemployment benefits
mmu = .15;
ttau = mmu*(1-aggEmployment)/aggEmployment;

% Aggregate Shocks
rrhoTFP = .859;										
ssigmaTFP = .014;

% Distribution of indv params log(lambda_i) ~ N(mu_l,-2*mu_l), 
% so lambda_i > 0 and E(lambda_i) = 1
mu_l = -.25; % Roughly calibrated to Piketty, Saez & Zucman (QJE 2018) Table I, log of P90-P20 ratio, post-tax

% Measurement error std. dev. of observed log output
ssigmaMeas = 0.02;


%% Set approximation parameters

global nEpsilon nAssets nAssetsFine nAssetsQuadrature ...
	nMeasure maxIterations tolerance dampening splineOpt displayOpt;

% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% MUST BE SET TO 0

% Whether to print out results from steady state computation
displayOpt = 'off';       % 'iter-detailed' or 'off'

% Order of approximation
nEpsilon = 2;
nAssets = 25; % number of gridpoints in spline approximation or polynomials in polynomial approximation

% Finer grid for analyzing policy functions
nAssetsFine = 100;

% Approximation of distribution
nMeasure = 3;
nAssetsQuadrature = 8;

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;


%% Save parameters

cd('./auxiliary_functions/dynare');
delete steady_vars.mat;
saveParameters;


%% Initial Dynare run

rng(rng_seed);
dynare firstOrderDynamics_polynomials noclearall nopathchange; % Run Dynare once to process model file


%% Simulate data

if is_data_gen == 0
    
    % Load previous data
    load('simul.mat')
    load('simul_data_hh_indv_param.mat');
    
else
    
    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw normalized individual incomes
    simul_data_hh = simulate_hh(sim_struct, ts_hh, N_hh);
    save('simul_data_hh.mat','simul_data_hh');
    
    % draw individual productivities and incomes
    simul_data_hh_indv_param = simulate_hh_indv_param(simul_data_hh);
    save('simul_data_hh_indv_param.mat','simul_data_hh_indv_param');
    
end


%% Compute likelihood

loglikes = nan(length(param1_vals),length(param2_vals));
loglikes_macro = nan(length(param1_vals),length(param2_vals));
loglikes_hh = nan(length(param1_vals),length(param2_vals));

disp('Computing likelihood...');
timer_likelihood = tic;

poolobj = parpool;

for iter_i=1:length(param1_vals) % For each parameter...
    
    for iter_j=1:length(param2_vals) % For each parameter...
        
        % Set new parameters
        bbeta = param1_vals(iter_i);
%         ssigmaMeas = param2_vals(iter_j);
        mu_l = param2_vals(iter_j);

%         fprintf(['%s' repmat('%6.4f ',1,2),'%s\n'], '[bbeta,ssigmaMeas] = [',bbeta,ssigmaMeas,']');
        fprintf(['%s' repmat('%6.4f ',1,2),'%s\n'], '[bbeta,mu_l] = [',bbeta,mu_l,']');

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state

        % Log likelihood of proposal
        [loglikes(iter_i,iter_j), loglikes_macro(iter_i,iter_j), loglikes_hh(iter_i,iter_j)] = ...
            loglike_compute_indv_param('simul.mat', simul_data_hh_indv_param, ts_hh, ...
                                       num_smooth_draws, num_interp, num_burnin_periods, ...
                                       M_, oo_, options_);
    
    end
    
end

delete(poolobj);

likelihood_elapsed = toc(timer_likelihood);
fprintf('%s%8.2f\n', 'Done. Elapsed minutes: ', likelihood_elapsed/60);

cd('../../');