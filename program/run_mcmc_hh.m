clear all;

model_name = 'hh';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));


%% Settings

% Decide what to do
is_run_dynare = true;   % Process Dynare model?
is_data_gen = true;     % Simulate data?
likelihood_type = 1;    % =1: macro + full-info micro; =2: macro only;
                        % =3: macro + 3 micro moments; =4: macro + 2 micro moments; =5: macro + 1 micro moment

% ID
serial_id = 1;          % ID number of current run (used in file names and RNG seeds)

% Model/data settings
T = 100;                % Number of periods of simulated macro data
ts_micro = 10:10:T;     % Time periods where we observe micro data
N_micro = 1e3;          % Number of households per non-missing time period

% File names
global mat_suff;
mat_suff = sprintf('%s%d%s%d%s%02d', '_N', N_micro, '_liktype', likelihood_type, '_', serial_id); % Suffix string for all saved .mat files
save_folder = fullfile(pwd, 'results'); % Folder for saving results

% Parameter transformation
if ismember(likelihood_type,[1 3 4]) % When mu_l is identified
    param_names = {'bbeta', 'ssigmaMeas', 'mu_l'};                      % Names of parameters to estimate
    transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2)) -exp(x(3))];     % Function mapping transformed parameters into parameters of interest
    param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2)) log(-x(3))];   % Function mapping parameters of interest into transformed parameters
else % When mu_l is not identified
    param_names = {'bbeta', 'ssigmaMeas'};                      % Names of parameters to estimate
    transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2))];     % Function mapping transformed parameters into parameters of interest
    param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2))];   % Function mapping parameters of interest into transformed parameters
end

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters

% Optimization settings
is_optimize = true;                             % Find posterior mode?
if ismember(likelihood_type,[1 3 4]) % When mu_l is identified
    [aux1, aux2, aux3] = meshgrid(linspace(0.8,0.99,5),linspace(0.001,0.05,5),linspace(-1,-0.01,5));
    optim_grid = [aux1(:), aux2(:), aux3(:)];   % Optimization grid
else % When mu_l is not identified
    [aux1, aux2] = meshgrid(linspace(0.8,0.99,5),linspace(0.001,0.05,5));
    optim_grid = [aux1(:), aux2(:)];
end
clearvars aux*;

% MCMC settings
if likelihood_type ~= 2
    mcmc_init = param_to_transf([.9 .06 -1]);   % Initial transformed draw (will be overwritten if is_optimize=true)
else % mu_l is not identified with macro data only
    mcmc_init = param_to_transf([.9 .06]);
end
mcmc_num_iter = 1e4;                  % Number of MCMC steps (total)
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
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20200813+serial_id;          % Random number generator seed
if likelihood_type == 1
    delete(gcp('nocreate'));    
    poolobj = parpool;                  % Parallel computing object
end

% Dynare settings
dynare_model = 'firstOrderDynamics_polynomials'; % Dynare model file


%% Calibrate parameters, execute initial Dynare processing

run_calib_dynare;


%% Simulate data

run_sim;


%% Measurement error

compute_meas_err_const; % Part of cov matrix of sample moments that doesn't change over parameter values


%% Find approximate mode

% Log likelihood function
ll_fct = @(M_, oo_, options_) aux_ll(simul_data_micro, ts_micro, ...
                                num_smooth_draws, num_burnin_periods, ...
                                num_interp, likelihood_type, ...
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

