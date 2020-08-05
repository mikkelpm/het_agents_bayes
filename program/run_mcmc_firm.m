clear all;
addpath(genpath('./functions'));
addpath(genpath('./firm_model/auxiliary_functions'));


%% Settings

% Decide what to do
is_run_dynare = false;  % Make Dynare process model files?
is_data_gen = true;     % Simulate data?
likelihood_type = 1;    % =1: Macro + full-info micro; =2: macro + full-info micro, no truncation; =3: macro + moments micro

% Model/data settings
T = 50;                                 % Number of periods of simulated macro data
ts_micro = 10:10:T;                     % Time periods where we observe micro data
N_micro = 1e2;                          % Number of micro entities per non-missing time period
trunc_logn = -1.1271;                   % Micro sample selection: Lower truncation point for log(n)

% Parameter transformation
param_names = {'rrhoProd', 'ssigmaProd'};                   % Names of parameters to estimate
transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2:end))];    % Function mapping transformed parameters into parameters of interest
param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2:end))];  % Function mapping parameters of interest into transformed parameters

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters

% Optimization settings
is_optimize = true;                                                    % Find posterior mode?
optim_grid = combvec(linspace(0.1,0.9,3),linspace(0.01,0.1,3))';    % Optimization grid

% MCMC settings
mcmc_init = param_to_transf([0.7 0.02]);% Distribution of initial transformed draw (will be overwritten if is_optimize=true)
mcmc_num_iter = 1e4;                    % Number of MCMC steps (total)
mcmc_thin = 1;                          % Store every X draws
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [50 200 500 1000];    % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = false;                % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                  % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)
mcmc_filename = 'firm';                 % File name of MCMC output

% Adaptive RWMH
mcmc_c = 0.55;
mcmc_ar_tg = 0.3;
mcmc_p_adapt = .95;

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 202006221;                   % Random number generator seed for initial simulation
poolobj = parpool;                      % Parallel computing object



%% Calibrate parameters and set numerical settings

run('firm_model/calibrate');

cd('./firm_model/dynare');
saveParameters;


%% Initial Dynare run

if is_run_dynare
    delete steady_vars.mat;
    dynare dynamicModel noclearall nopathchange; % Run Dynare once to process model file
    save dynare.mat M_ oo_ options_;
else
    load('dynare.mat');
end


%% Simulate data

set_dynare_seed(rng_seed);  % Seed Dynare RNG
rng(rng_seed, 'twister');   % Seed Matlab RNG

if ~is_data_gen
    % Load previous data
    load('simul_data_micro.mat');
else
    % Simulate
    simul_data;
end


%% Find approximate mode

% Log likelihood function
ll_fct = @(M_, oo_, options_) aux_ll(simul_data_micro, ts_micro, ...
                                num_smooth_draws, num_burnin_periods, ...
                                trunc_logn, likelihood_type, ...
                                M_, oo_, options_);

% Optimization
if is_optimize
    approx_mode;
end


%% Run MCMC iterations

mcmc_iter;


%% Save results

cd('../');
addpath('results');
save(fullpath('results', mcmc_filename));

delete(poolobj);


%% Auxiliary function

function [the_loglike, the_loglike_macro, the_loglike_micro] = ...
         aux_ll(simul_data_micro, ts_micro, ...
                num_smooth_draws, num_burnin_periods, ...
                trunc_logn, likelihood_type, ...
                M_, oo_, options_)

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
        compute_meas_err;       % Update measurement error var-cov matrix for sample moments

        % Log likelihood of proposal
        switch likelihood_type
            case 1 % Macro + full info micro
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                                   num_smooth_draws, trunc_logn, num_burnin_periods, ...
                                   M_, oo_, options_);
            case 2 % Macro + full info micro, ignore truncation
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                                   num_smooth_draws, -Inf, num_burnin_periods, ...
                                   M_, oo_, options_);
            case 3 % Macro + moments w/ SS meas. err.
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute('simul_moments.mat', [], ts_micro, ...
                                   0, [], num_burnin_periods, ...
                                   M_, oo_, options_);
        end
        
end