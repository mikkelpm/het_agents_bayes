clear all;
addpath(genpath('./functions'));

model_name = 'hh';
addpath(genpath(['./' model_name '_model/auxiliary_functions']));


%% Settings

% Decide what to do
is_run_dynare = true;   % Process Dynare model?
is_data_gen = true;     % Simulate data?
likelihood_type = 1;    % =1: Macro + full-info micro; =2: macro + full-info micro, no truncation; =3: macro + moments micro; =4: macro only

% Model/data settings
T = 100;                % Number of periods of simulated macro data
ts_micro = 10:10:T;     % Time periods where we observe micro data
N_micro = 1e2;          % Number of households per non-missing time period

% Parameter transformation
param_names = {'bbeta', 'ssigmaMeas', 'mu_l'};                      % Names of parameters to estimate
transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2)) -exp(x(3))];     % Function mapping transformed parameters into parameters of interest
param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2)) log(-x(3))];   % Function mapping parameters of interest into transformed parameters

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters

% Optimization settings
is_optimize = true;                         % Find posterior mode?
[aux1, aux2, aux3] = meshgrid(linspace(0.8,0.99,3),linspace(0.001,0.05,3),linspace(-1,-0.01,3));
optim_grid = [aux1(:), aux2(:), aux3(:)];   % Optimization grid

% MCMC settings
mcmc_init = param_to_transf([.9 .06 -1]);   % Initial transformed draw (will be overwritten if is_optimize=true)
mcmc_num_iter = 1e4;                  % Number of MCMC steps (total)
mcmc_thin = 1;                          % Store every X draws
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [50 200 500 1000];    % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = false;                % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                  % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)
mcmc_filename = [model_name '_liktype' num2str(likelihood_type) '_N' num2str(N_micro) '_']; % File name of MCMC output

% for adaptive RWMH
mcmc_c = 0.55;
mcmc_ar_tg = 0.3;
mcmc_p_adapt = .95;

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 202006221;                   % Random number generator seed for initial simulation
poolobj = parpool;                      % Parallel computing object

global mat_suff;
mat_suff = sprintf('%02d', 1);          % Suffix string for all saved .mat files


%% Calibrate parameters and set numerical settings

run([model_name '_model/calibrate']);

cd(['./' model_name '_model/dynare']);
saveParameters;


%% Initial Dynare processing

if is_run_dynare
    dynare firstOrderDynamics_polynomials noclearall nopathchange; % Run Dynare once to process model file
else
    load('firstOrderDynamics_polynomials_results');
    check_matlab_path(false);
    dynareroot = dynare_config(); % Add Dynare sub-folders to path
end


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
                                num_interp, likelihood_type, ...
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
save_mat(fullfile('results', mcmc_filename));

delete(poolobj);


%% Auxiliary function

function [the_loglike, the_loglike_macro, the_loglike_micro] = ...
         aux_ll(data_micro, ts_micro, ...
                num_smooth_draws, num_burnin_periods, ...
                num_interp, likelihood_type, ...
                M_, oo_, options_)
        
        global mat_suff;

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
        compute_meas_err;       % Update measurement error var-cov matrix for sample moments
        
        % Micro likelihood function
        aaBar_local = aaBar;
        mmu_local = mmu;
        ttau_local = ttau;
        mu_l_local = mu_l;
        num_mom = 3;
        likelihood_micro_fct = @(smooth_draw,it) ...
                               likelihood_micro(smooth_draw, ts_micro(it), data_micro, it, ...
                                                aaBar_local, mmu_local, ttau_local, mu_l_local, num_mom, ...
                                                num_interp);
        
        % Macro state variables used in micro likelihood
        smooth_vars = [{'w'; 'r'; 'lag_mHat_1' ; 'lag_mHat_2'};
                       str_add_numbers('lag_moment_1_', 1:num_mom);
                       str_add_numbers('lag_moment_2_', 1:num_mom);
                       str_add_numbers('measureCoefficient_1_', 1:num_mom);
                       str_add_numbers('measureCoefficient_2_', 1:num_mom)];
        
        % Log likelihood computation
        switch likelihood_type
            case 1 % Macro + full info micro
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, num_smooth_draws, ...
                                   M_, oo_, options_, ...
                                   ts_micro, @(smooth_draw,it) likelihood_micro_fct(smooth_draw,it));
            case 2 % Macro + up to 3rd moments w/ SS meas. err.
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul_moments', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], []);
            case 3 % Macro + up to 2nd moments w/ SS meas. err.
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul_moments2', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], []);
            case 4 % Macro only
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], []);
        end
        
end

function str_cell = str_add_numbers(str, numbers)
    str_cell = cellfun(@(x) sprintf('%s%d', str, x), num2cell(numbers), 'UniformOutput', false)';
end
