clear all;

model_name = 'firm';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));


%% Settings

% Decide what to do
is_run_dynare = true;   % Process Dynare model?
is_data_gen = true;     % Simulate data?
likelihood_type = 1;    % =1: macro + full-info micro; =2: macro + full-info micro, ignore truncation; =3: macro + micro moments

% ID
serial_id = 1;          % ID number of current run (used in file names and RNG seeds)

% Model/data settings
T = 50;                 % Number of periods of simulated macro data
ts_micro = 10:10:T;     % Time periods where we observe micro data
N_micro = 1e3;          % Number of micro entities per non-missing time period
trunc_quant = 0.9;      % Micro sample selection: Lower truncation quantile for labor (steady state distribution)

% Suffix string for all saved .mat files
global mat_suff;
mat_suff = sprintf('%s%d%s%d%s%02d', '_trunc', 100*trunc_quant, '_liktype', likelihood_type, '_', serial_id);

% Parameter transformation
param_names = {'rrhoProd', 'ssigmaProd'};               % Names of parameters to estimate
n_param = length(param_names);
transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2))];    % Function mapping transformed parameters into parameters of interest
param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2))];  % Function mapping parameters of interest into transformed parameters

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters

% Optimization settings
is_optimize = true;                 % Find posterior mode?
[aux1, aux2] = meshgrid(linspace(0.1,0.9,5),linspace(0.01,0.1,5));
optim_grid = [aux1(:), aux2(:)];    % Optimization grid
clearvars aux1 aux2;

% MCMC settings
mcmc_init = param_to_transf([0.7 0.02]);% Initial transformed draw (will be overwritten if is_optimize=true)
mcmc_num_iter = 1e4;                    % Number of MCMC steps (total)
mcmc_thin = 1;                          % Store every X draws
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [50 200 500 1000];    % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = false;                % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                  % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)

% Adaptive RWMH
mcmc_c = 0.55;                          % Updating rate parameter
mcmc_ar_tg = 0.3;                       % Target acceptance rate
mcmc_p_adapt = .95;                     % Probability of non-diffuse proposal

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20200813+serial_id;          % Random number generator seed
if likelihood_type == 1
    delete(gcp('nocreate'));    
    poolobj = parpool;                  % Parallel computing object
end


%% Calibrate parameters and set numerical settings

run([model_name '_model/calibrate']);

cd(['./' model_name '_model/dynare']);
saveParameters;


%% Initial Dynare processing

if is_run_dynare
    dynare dynamicModel noclearall nopathchange; % Run Dynare once to process model file
else
    load('dynamicModel_results');
    check_matlab_path(false);
    dynareroot = dynare_config(); % Add Dynare sub-folders to path
end


%% Truncation point

compute_steady_state; % Re-compute steady state
trunc_logn = M_.steady_vars.smpl_m1+norminv(trunc_quant)*sqrt(M_.steady_vars.smpl_m3); % Lower truncation value for log(n)


%% Simulate data

set_dynare_seed(rng_seed);  % Seed Dynare RNG
rng(rng_seed, 'twister');   % Seed Matlab RNG

if ~is_data_gen
    % Load previous data
    load_mat('simul_data_micro');
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

cd('../../');
mkdir('results');
save_mat(fullfile('results', model_name));

if likelihood_type == 1
    delete(gcp('nocreate'));
end

%% Auxiliary likelihood function

function [the_loglike, the_loglike_macro, the_loglike_micro] = ...
         aux_ll(data_micro, ts_micro, ...
                num_smooth_draws, num_burnin_periods, ...
                trunc_logn, likelihood_type, ...
                M_, oo_, options_)
        
        global mat_suff;

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
        compute_meas_err;       % Update measurement error var-cov matrix for sample moments
        
        % Macro state variables used in micro likelihood
        smooth_vars = [{'logWage'; 'aggregateTFP'};
                       str_add_numbers('lag_moment_', 1:5)];
        
        % Parameters passed to micro likelihood function
        param = [nnu ttheta trunc_logn];
        
        % Log likelihood computation
        switch likelihood_type
            case 1 % Macro + full info micro
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, num_smooth_draws, ...
                                   M_, oo_, options_, ...
                                   data_micro, ts_micro, param);
            case 2 % Macro + full info micro, ignore truncation
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, num_smooth_draws, ...
                                   M_, oo_, options_, ...
                                   data_micro, ts_micro, [param(1:2) -Inf]);
            case 3 % Macro + moments w/ SS meas. err.
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul_moments', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], [], param);
        end
        
end