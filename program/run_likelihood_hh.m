clear all;

model_name = 'hh';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));


%% Settings

% Decide what to do
is_run_dynare = false;   % Process Dynare model?
is_data_gen = true;     % Simulate data?

% ID
serial_id = 1;          % ID number of current run (used in file names and RNG seeds)

% Model/data settings
T = 100;                % Number of periods of simulated macro data
ts_micro = 10:10:T;     % Time periods where we observe micro data
N_micro = 1e3;          % Number of households per non-missing time period

% File names
global mat_suff;
mat_suff = sprintf('%s%d%s%02d', '_likelihood_N', N_micro, '_', serial_id); % Suffix string for all saved .mat files
save_folder = fullfile(pwd, 'results'); % Folder for saving results

% Parameters to evaluate
param_names = {'bbeta', 'ssigmaMeas', 'mu_l'};          % Names of parameters to evaluate
param_vals_mult = [0.75 1.25];                          % Lowest and highest multiples of true value on fine grid
param_vals_num_fine = 100;                              % Number of fine grid points (excluding true values) in interval defined by "param_vals_mult"
param_space =        [0   1;   0   Inf; -Inf 0];        % Boundaries of parameter space to enforce on fine grid
param_vals_include = [nan nan; nan 0.3; -0.5 nan];      % Include these lower/upper endpoints for each parameter (nan: no override)
param_vals_num_coarse = 50;                             % Number of coarse grid points outside interval defined by "param_vals_mult"

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20200813+serial_id;          % Random number generator seed

delete(gcp('nocreate'));    
poolobj = parpool;                      % Parallel computing object


%% Calibrate parameters and set numerical settings

run([model_name '_model/calibrate']);

cd(['./' model_name '_model/dynare']);
saveParameters;
economicParameters_true = load_mat('economicParameters'); % Store true parameters


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


%% Parameter combinations to evaluate

n_param = length(param_names);
params_truth = nan(1,n_param);

% True parameter values
for i_param = 1:n_param
    params_truth(i_param) = eval(param_names{i_param});
end

% Combinations of parameters for likelihood evaluation
% Vary each parameter individually, keeping other parameters at true values
lik_grid = nan(0,n_param);
len_lik = nan(1,n_param); % Length of lik_grid for each parameter
for i_param = 1:n_param

    % Fine grid for parameter
    aux = params_truth(i_param)*param_vals_mult;
    the_grid = sort([params_truth(i_param) ...
                     linspace(max(min(aux),param_space(i_param,1)), ...
                              min(max(aux),param_space(i_param,2)), ...
                              param_vals_num_fine) ...
                    ]);
    
    % Extend with coarse grid to left/right, if desired
    if param_vals_include(i_param,1)<the_grid(1)
        the_lin = linspace(param_vals_include(i_param,1), the_grid(1), param_vals_num_coarse+1);
        the_grid = [the_lin(1:end-1) the_grid];
    end
    if param_vals_include(i_param,2)>the_grid(end)
        the_lin = linspace(the_grid(end), param_vals_include(i_param,2), param_vals_num_coarse+1);
        the_grid = [the_grid the_lin(2:end)];
    end
    
    % Add to grid of all parameter combinations
    len_lik(i_param) = length(the_grid);
    aux2 = repmat(params_truth, len_lik(i_param), 1);
    aux2(:,i_param) = the_grid;
    lik_grid = [lik_grid; aux2];
    
end
clearvars aux aux2 the_grid the_lin;


%% Evaluate likelihoods

lik_numgrid = size(lik_grid,1);
lik_all = nan(lik_numgrid,3,5);

fprintf('\nLikelihoood...\n');
timer_lik = tic;

for i_lik=1:lik_numgrid % Cycle through parameters

    the_param = lik_grid(i_lik,:);
    print_param(the_param, param_names, 'current');
    update_param(the_param, param_names);
    
    saveParameters;         % Save parameter values to files
    setDynareParameters;    % Update Dynare parameters in model struct
    try
        compute_steady_state;   % Compute steady state
    catch ME
        disp('Error encountered in steady state computation. Message:');
        disp(ME.message);
        continue;
    end

    for i_type=1:5 % Cycle through likelihood types
        
        fprintf('%s%d\n', 'Likelihood type: ', i_type);
        
        try
            [lik_all(i_lik,1,i_type),lik_all(i_lik,2,i_type),lik_all(i_lik,3,i_type)] = ...
                                    aux_ll(simul_data_micro, ts_micro, ...
                                    num_smooth_draws, num_burnin_periods, ...
                                    num_interp, i_type, ...
                                    M_, oo_, options_); % Only compute steady state once
        catch ME
            disp('Error encountered in likelihoood computation. Message:');
            disp(ME.message);
        end
        
    end

    % Print progress
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_lik, '/', lik_numgrid);

end

%% Save results

mkdir(save_folder);
save_mat(fullfile(save_folder, model_name));

delete(gcp('nocreate'));


%% Auxiliary function

function [the_loglike, the_loglike_macro, the_loglike_micro] = ...
         aux_ll(data_micro, ts_micro, ...
                num_smooth_draws, num_burnin_periods, ...
                num_interp, likelihood_type, ...
                M_, oo_, options_)
        
        global mat_suff;
        
        saveParameters;         % Save parameter values to files
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
