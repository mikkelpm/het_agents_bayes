clear all;

% Simulate and estimate heterogeneous firm model

model_name = 'firm';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));


%% Settings

% Decide what to do
is_run_dynare = true;   % Process Dynare model?
is_data_gen = true;     % Simulate data?
likelihood_type = 1;    % =1: macro + full-info micro; =2: macro only

% Experiment type
exper_type = 1; % =1: estimate adjustment cost parameters, no truncation
                % =2: estimate adjustment cost parameters, with truncation
                % =3: estimate idiosyncratic productivity parameters, no truncation

% ID
serial_id = 1;          % ID number of current run (used in file names and RNG seeds)

% Model/data settings
T = 50;                 % Number of periods of simulated macro data
if exper_type < 3
    ts_micro = 1:T;     % Time periods when we observe micro data
else
    ts_micro = 10:10:T;
end
if exper_type ~= 2
    N_micro = 1e3;      % Number of micro entities per non-missing time period
    trunc_quant = 0;  	% Micro sample selection: Lower truncation quantile for labor (steady state distribution), =0 if no truncation
else
    N_micro = 1e4;
    trunc_quant = 0.9;
end

% File names
global mat_suff;
mat_suff = sprintf('%s%d%s%d%s%02d', '_exp', exper_type, '_liktype', likelihood_type, '_', serial_id); % Suffix string for all saved .mat files
save_folder = fullfile(pwd, 'results'); % Folder for saving results

% Parameter transformations and prior
if exper_type < 3
    param_names = {'ppsiCapital', 'aaUpper'};   % Names of parameters to estimate (NOTE: the function "setParameters" sets aaLower=-aaUpper)
    transf_to_param = @(x) exp(x);              % Function mapping transformed parameters into parameters of interest
    param_to_transf = @(x) log(x);              % Function mapping parameters of interest into transformed parameters
    prior_logdens_transf = @(x) sum(x);         % Log prior density of transformed parameters
else
    param_names = {'rrhoProd', 'ssigmaProd'};
    transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2))];
    param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2))];
    prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));
end

% Optimization settings
is_optimize = true;                 % Find posterior mode?
if exper_type < 3
    [aux1, aux2] = meshgrid(linspace(0.001,0.05,20),linspace(0.001,0.05,20));
else
    [aux1, aux2] = meshgrid(linspace(0.1,0.9,5),linspace(0.01,0.1,5));
end
optim_grid = [aux1(:), aux2(:)];    % Optimization grid
clearvars aux1 aux2;

% MCMC settings
if exper_type < 3
    mcmc_init = param_to_transf([0.01 0.005]); % Initial transformed draw (will be overwritten if is_optimize=true)
else
    mcmc_init = param_to_transf([0.7 0.02]);
end
mcmc_num_iter = 1e4;                    % Number of MCMC steps (total)
mcmc_thin = 1;                          % Store every X draws
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [50 200 500 1000];    % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = false;                % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                  % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)

% Adaptive RWMH
mcmc_c = 0.55;                          % Updating rate parameter
mcmc_ar_tg = 0.3;                       % Target acceptance rate
if exper_type < 3
    mcmc_p_adapt = 1;                   % Probability of non-diffuse proposal
else
    mcmc_p_adapt = .95;
end

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20200813+serial_id;          % Random number generator seed
if likelihood_type == 1
    delete(gcp('nocreate'));
    poolobj = parpool;                  % Parallel computing object
end

% Dynare settings
dynare_model = 'dynamicModel';          % Dynare model file


%% Calibrate parameters, execute initial Dynare processing

run_calib_dynare;
% Approximates cross-sec distribution of firm productivity and capital
% using multivariate normal distribution, as in Winberry (2018)


%% Truncation point

compute_steady_state; % Re-compute steady state
empl_ss_mean = -(M_.steady_vars.moment_1-log(M_.steady_vars.wage)+log(nnu)+ttheta*M_.steady_vars.moment_2)/(nnu-1); % Steady-state cross-sectional mean of log employment
empl_ss_var = (M_.steady_vars.moment_3+ttheta^2*M_.steady_vars.moment_5+2*ttheta*M_.steady_vars.moment_4)/(nnu-1)^2; % Steady-state cross-sectional variance of log employment
trunc_logn = empl_ss_mean+norminv(trunc_quant)*sqrt(empl_ss_var); % Lower truncation value for log(n), exploiting approximating normal distribution


%% Simulate data

run_sim;


%% Find approximate mode

% Log likelihood function
ll_fct = @(M_, oo_, options_) aux_ll(simul_data_micro, ts_micro, ...
                                num_smooth_draws, num_burnin_periods, ...
                                trunc_logn, likelihood_type, ...
                                M_, oo_, options_, ...
                                true);

% Optimization
if is_optimize
    approx_mode;
end


%% Run MCMC iterations

mkdir(save_folder);
mcmc_iter;


%% Save results

save_mat(fullfile(save_folder, model_name));

if likelihood_type == 1
    delete(poolobj);
end

