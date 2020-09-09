% based on run_memc_hh_array.m & approx_mode.m

clear all;

model_name = 'hh';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));

array_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));


%% Settings

% Decide what to do
is_run_dynare = false;   % Process Dynare model?
is_data_gen = false;     % Simulate data?
likelihood_type = array_id; % floor((array_id-1)/10);    % =1: macro + full-info micro; =2: macro only;
                        % =3: macro + 3 micro moments; =4: macro + 2 micro moments; =5: macro + 1 micro moment

% ID
serial_id = 1; % array_id-likelihood_type*10;          % ID number of current run (used in file names and RNG seeds)

% Model/data settings
T = 100;                % Number of periods of simulated macro data
ts_micro = 10:10:T;     % Time periods where we observe micro data
N_micro = 1e3;          % Number of households per non-missing time period

% Suffix string for all saved .mat files
global mat_suff;
mat_suff = sprintf('%s%d%s%d%s%02d', '_N', N_micro, '_liktype', likelihood_type, '_', serial_id);

% Parameter transformation
if likelihood_type ~= 2
    param_names = {'bbeta', 'ssigmaMeas', 'mu_l'};                      % Names of parameters to estimate
    n_param = length(param_names);
    transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2)) -exp(x(3))];     % Function mapping transformed parameters into parameters of interest
    param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2)) log(-x(3))];   % Function mapping parameters of interest into transformed parameters
else % mu_l is not identified with macro data only
    param_names = {'bbeta', 'ssigmaMeas'};                      % Names of parameters to estimate
    n_param = length(param_names);
    transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2))];     % Function mapping transformed parameters into parameters of interest
    param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2))];   % Function mapping parameters of interest into transformed parameters
end

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters

% Likelihood parameter settings
param_vals_mult = unique([1 linspace(0.75,1.25,51)]); % Multiples of true parameter to compute
n_lik = length(param_vals_mult);
params_truth = [0.96 0.02 -0.25];
params_truth = params_truth(:,1:n_param);   % Likelihood parameters
lik_grid = repmat(params_truth,n_lik*n_param,1);   % Likelihood parameters
for i_param = 1:n_param
    lik_grid((i_param-1)*n_lik+(1:n_lik),i_param) = params_truth(i_param)*param_vals_mult;
end

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


%% Evaluate likelihoods

% Evaludate likelihoods

lik_numgrid = size(lik_grid,1);
lik_all = nan(lik_numgrid,3);

fprintf('\nLikelihoood...\n');
timer_lik = tic;

for i_lik=1:lik_numgrid % Cycle through parameters

    the_param = lik_grid(i_lik,:);
    print_param(the_param, param_names, 'current');
    for i=1:n_param
        eval(sprintf('%s%s%f%s', param_names{i}, '=', the_param(i), ';'));
    end

    try
        [lik_all(i_lik,1),lik_all(i_lik,2),lik_all(i_lik,3)] = ...
                                aux_ll(simul_data_micro, ts_micro, ...
                                num_smooth_draws, num_burnin_periods, ...
                                num_interp, likelihood_type, ...
                                M_, oo_, options_);
    catch ME
        disp('Error encountered. Message:');
        disp(ME.message);
    end

    % Print progress
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_lik, '/', lik_numgrid);

end

%% Save results

cd('../../');
mkdir('results');
save_mat(fullfile('results', [model_name '_likelihood']));

if likelihood_type == 1
    delete(gcp('nocreate'));
end


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
        
        % Macro state variables used in micro likelihood
        num_mom = 3;
        smooth_vars = [{'w'; 'r'; 'lag_mHat_1' ; 'lag_mHat_2'};
                       str_add_numbers('lag_moment_1_', 1:num_mom);
                       str_add_numbers('lag_moment_2_', 1:num_mom);
                       str_add_numbers('measureCoefficient_1_', 1:num_mom);
                       str_add_numbers('measureCoefficient_2_', 1:num_mom)];
        
        % Parameters passed to micro likelihood function
        param = [aaBar mmu ttau mu_l num_mom num_interp];
        
        % Log likelihood computation
        switch likelihood_type
            case 1 % Macro + full info micro
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, num_smooth_draws, ...
                                   M_, oo_, options_, ...
                                   data_micro, ts_micro, param);
            case 2 % Macro only
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], [], param);
            case 3 % Macro + up to 3rd moments
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul_moments', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], [], param);
            case 4 % Macro + up to 2nd moments
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul_moments2', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], [], param);
            case 5 % Macro + only 1st moments
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul_moments1', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], [], param);
        end
        
end
