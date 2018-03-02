clear all;

%% Settings

is_data_gen = 1; % whether simulate data
is_profile = 0; %whether run profiler for execution time

bbetas = linspace(0.9,0.99,10);     % beta values to loop over

T = 200;                            % Number of periods of simulated macro data
ts_hh = 20:20:200;                  % Time periods where we observe micro data
N_hh = 1e3;                         % Number of households per non-missing time period

constr_tol = 1e-6;                  % Numerical tolerance for whether assets are at borrowing constraint
num_burnin_periods = 100;           % Number of burn-in periods for simulations
num_smooth_draws = 25;              % Number of draws from the smoothing distribution (for each beta)

rng_seed = 201803021;                % Random number generator seed for initial simulation

tag_date = datestr(now,'yyyymmdd');


%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP;
	
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
saveParameters;


%% Initial Dynare run

dynare firstOrderDynamics_polynomials noclearall;                   % Run Dynare once to process model file


%% Simulate data

if is_data_gen
    
    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw individual incomes
    simul_data_hh = simulate_hh(sim_struct, ts_hh, N_hh);
    save('simul_data_hh.mat','simul_data_hh');
    
else
    
    load('simul.mat')
    load('simul_data_hh.mat');
    
end


%% Smoothing and likelihood

% Loop over beta values
loglikes = zeros(size(bbetas));
loglikes_macro = zeros(size(bbetas));
loglikes_hh = zeros(size(bbetas));

for i_beta=1:length(bbetas) % For each alpha...

    fprintf('%s%6.4f\n', 'beta=', bbetas(i_beta));
    bbeta = bbetas(i_beta);                          % Set beta
    saveParameters;
    setDynareParameters;

    [loglikes(i_beta), loglikes_macro(i_beta), loglikes_hh(i_beta)] = ...
        loglike_compute('simul.mat', simul_data_hh, ts_hh, num_smooth_draws, num_burnin_periods, constr_tol, M_, oo_, options_);
    
end

cd('../');

if is_profile
    profsave(profile('info'),['profile_results_' tag_date])
end