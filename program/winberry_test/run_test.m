clear all;

%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data
is_profile = 0; %whether run profiler for execution time

% Model/data settings
T = 100;                            % Number of periods of simulated macro data
ts_hh = 20:20:100;                  % Time periods where we observe micro data
N_hh = 1e2;                         % Number of households per non-missing time period

% Prior
prior_logdens_logbeta = @(x) x;     % Log prior density of log(beta)
prior_init_logbeta = @() log(0.9);     % Distribution of initial log(beta) draw

% MCMC settings
mcmc_num_draws = 1000;              % Number of MCMC steps (total, including initial phase)
mcmc_num_adapt = 50;                % Number of initial steps with larger (but gradually decreasing) step size
mcmc_stepsize = 1e-3;               % MCMC step size after initial phase
mcmc_stepsize_init = 10*mcmc_stepsize; % Initial MCMC steps size (gradually decreased during initial phase)
mcmc_smooth_draws = 500;            % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
mcmc_filename = 'output/mcmc.mat';  % File name of MCMC output

% Numerical settings
constr_tol = 1e-6;                  % Numerical tolerance for whether assets are at borrowing constraint
num_burnin_periods = 100;           % Number of burn-in periods for simulations
rng_seed = 20180305;                % Random number generator seed for initial simulation

% Profiler save settings
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

rng(rng_seed);
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


%% MCMC

curr_logbeta = prior_init_logbeta();        % Initial draw
post_draws_beta = nan(mcmc_num_draws,1);
accepts = zeros(mcmc_num_draws,1);

curr_loglike = -Inf;
loglikes_prop = nan(mcmc_num_draws,1);
loglikes_prop_macro = nan(mcmc_num_draws,1);
loglikes_prop_hh = nan(mcmc_num_draws,1);

disp('MCMC...');
timer_mcmc = tic;

for i_mcmc=1:mcmc_num_draws % For each MCMC step...

    fprintf('%s%6.4f\n', 'current  beta=', exp(curr_logbeta));
    
    % Proposed log(beta) (modified to always start with initial draw)
    the_stepsize = (i_mcmc>1)*(mcmc_stepsize + max(1-(i_mcmc-2)/mcmc_num_adapt,0)*(mcmc_stepsize_init-mcmc_stepsize)); % Step size
    prop_logbeta = curr_logbeta + the_stepsize*randn(); % Proposal
    fprintf('%s%6.4f\n', 'proposed beta=', exp(prop_logbeta));
    
    % Set new beta
    bbeta = exp(prop_logbeta);
    saveParameters;
    setDynareParameters;

    % Log likelihood of proposal
    success = true;
    try
        [loglikes_prop(i_mcmc), loglikes_prop_macro(i_mcmc), loglikes_prop_hh(i_mcmc)] = ...
            loglike_compute('simul.mat', simul_data_hh, ts_hh, mcmc_smooth_draws, num_burnin_periods, constr_tol, M_, oo_, options_);
    catch ME
        success = false;
        disp('Error encountered in likelihood evaluation. Message:');
        disp(ME.message);
    end
    
    if success
        % Log prior ratio
        logratio_prior = prior_logdens_logbeta(prop_logbeta)-prior_logdens_logbeta(curr_logbeta);

        % Accept/reject
        if loglikes_prop(i_mcmc)-curr_loglike+logratio_prior > log(rand())
            fprintf('%s\n', 'Accepted!');
            curr_logbeta = prop_logbeta;
            curr_loglike = loglikes_prop(i_mcmc);
            accepts(i_mcmc) = 1;
        end
    end
    
    % Store
    post_draws_beta(i_mcmc) = exp(curr_logbeta);
    
    % Print acceptance rate
    fprintf('%s%5.1f%s\n', 'Accept. rate last 100: ', 100*mean(accepts(max(i_mcmc-99,1):i_mcmc)), '%');
    fprintf('%s%5.1f%s\n\n', 'Progress: ', 100*i_mcmc/mcmc_num_draws, '%');
    
end

mcmc_elapsed = toc(timer_mcmc);
fprintf('%s%8.2f\n', 'MCMC done. Elapsed minutes: ', mcmc_elapsed/60);

save(mcmc_filename, 'mcmc_*', 'post_draws_beta', 'accepts', 'loglikes_prop_*');

cd('../');

if is_profile
    profsave(profile('info'),['profile_results_' tag_date])
end