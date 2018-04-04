clear all;

%% Settings

is_data_gen = 0; % whether simulate data:  
                 % 0: no simulation; 
                 % 1: simulation (with indv params, based on simulation without indv params)
                 % 2: simulation (with indv params, start from scratch)
is_profile = 0; % whether run profiler for execution time

bbetas = linspace(0.93,0.99,3);     % beta values to loop over
mu_ls = [1 2 3];           % params for dist of indv params

T = 200;                            % Number of periods of simulated macro data
ts_hh = 20:20:200;                  % Time periods where we observe micro data
N_hh = 1e2;                         % Number of households per non-missing time period

num_burnin_periods = 100;           % Number of burn-in periods for simulations
num_smooth_draws = 500;              % Number of draws from the smoothing distribution (for each beta)

rng_seed = 20180329;                % Random number generator seed for initial simulation

tag_date = datestr(now,'yyyymmdd');


%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu ttau rrhoTFP ssigmaTFP mu_l;
	
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

% Distribution of indv params log(lambda_i) ~ N(-1/2,1), 
% so lambda_i > 0 and E(lambda_i) = 1
mu_l = 2;

%% Set approximation parameters

global nEpsilon nAssets nAssetsFine nAssetsQuadrature ...
	nMeasure maxIterations tolerance dampening splineOpt displayOpt;

% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% if splineOpt = 1, use splines to approximate savings policy; if splineOpt = 0, use polynomials
				% to approximate conditional expectation function

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

cd('./Auxiliary Functions');
delete steady_vars.mat;
saveParameters;


%% Initial Dynare run

dynare firstOrderDynamics_polynomials noclearall;                   % Run Dynare once to process model file


%% Simulate data

if is_data_gen == 0
    
    load('simul.mat')
    load('simul_data_hh_indv_param.mat');
    
elseif is_data_gen == 1
    
    load('simul.mat')
    load('simul_data_hh.mat');
    
    % draw individual productivities and incomes
    simul_data_hh_indv_param = simulate_hh_indv_param(simul_data_hh);
    save('simul_data_hh_indv_param.mat','simul_data_hh_indv_param');
    
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


%% Smoothing and likelihood

n_beta = length(bbetas);
n_mul = length(mu_ls);

% Loop over beta values
loglikes = zeros(n_mul,n_beta);
loglikes_macro = zeros(n_mul,n_beta);
loglikes_hh = zeros(n_mul,n_beta);

poolobj = parpool;

for i_beta=1:n_beta % For each param...
    
    bbeta = bbetas(i_beta); % Set beta
    saveParameters;         % Save parameter values to files
    setDynareParameters;    % Update Dynare parameters in model struct
    compute_steady_state;   % Compute steady state once and for all
    
    for i_mul = 1:n_mul
        mu_l = mu_ls(i_mul);
        fprintf([repmat('%s%6.4f',1,2) '\n'], 'beta=', bbeta, ', mu_l=', mu_l);
        
        [loglikes(i_mul,i_beta), loglikes_macro(i_mul,i_beta), loglikes_hh(i_mul,i_beta)] = ...
            loglike_compute_indv_param('simul.mat', simul_data_hh_indv_param, ts_hh, num_smooth_draws, num_burnin_periods, M_, oo_, options_);
    end
    
end

delete(poolobj);

cd('../');

if is_profile
    profsave(profile('info'),['profile_results_' tag_date]);
end