clear all;

% Compute likelihood function (for various observables) in heterogeneous household model

model_name = 'firm';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));


%% Settings

% Decide what to do
is_run_dynare = false;   % Process Dynare model?
is_data_gen = true;     % Simulate data?

% ID
serial_id = 1;          % ID number of current run (used in file names and RNG seeds)

% Model/data settings
T = 50;                % Number of periods of simulated macro data
ts_micro = 1:T;     % Time periods where we observe micro data
N_micro = 1e3;          % Number of households per non-missing time period

% File names
global mat_suff;
mat_suff = sprintf('%s%02d', '_likelihood_', serial_id); % Suffix string for all saved .mat files
save_folder = fullfile(pwd, 'results2'); % Folder for saving results

% Parameters to evaluate
param_names = {'ppsiCapital', 'aaUpper'};          % Names of parameters to evaluate
param_vals_mult = [0.5 1.5];                          % Lowest and highest multiples of true value on fine grid
param_vals_num = 20;                              % Number of fine grid points (excluding true values) in interval defined by "param_vals_mult"
param_space =        [0   Inf;   0   Inf];        % Boundaries of parameter space to enforce on fine grid

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20200813+serial_id;          % Random number generator seed
delete(gcp('nocreate'));    
poolobj = parpool;                      % Parallel computing object

% Dynare settings
dynare_model = 'dynamicModel'; % Dynare model file


%% Calibrate parameters, execute initial Dynare processing

run_calib_dynare;


%% Simulate data

trunc_logn = -Inf;
run_sim;


%% Parameter combinations to evaluate

n_param = length(param_names);
params_truth = nan(1,n_param);

% True parameter values
for i_param = 1:n_param
    params_truth(i_param) = eval(param_names{i_param});
end

param_vals = linspace(param_vals_mult(1), param_vals_mult(2), param_vals_num)' * params_truth;

% Combinations of parameters for likelihood evaluation
lik_grid = [repmat(param_vals(:,1), param_vals_num, 1) kron(param_vals(:,2), ones(param_vals_num,1))]; % Parameter value combinaitions on which the likelihoods are computed


%% Evaluate likelihoods

% Simulate random shocks that are held fixed across parameters
sim_shocks = simulate_shocks(M_, T+num_burnin_periods, num_smooth_draws);

% Likelihood

lik_numgrid = size(lik_grid,1);
lik_all = nan(lik_numgrid,3,1);

fprintf('\nLikelihood...\n');
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

    for i_type=1:1 % Cycle through likelihood types
        
        fprintf('%s%d\n', 'Likelihood type: ', i_type);
        
        try
            [lik_all(i_lik,1,i_type),lik_all(i_lik,2,i_type),lik_all(i_lik,3,i_type)] = ...
                         aux_ll(simul_data_micro, ts_micro, ...
                                num_smooth_draws, num_burnin_periods, ...
                                -Inf, i_type, ...
                                M_, oo_, options_, ...
                                false, ...   % Only compute steady state once
                                sim_shocks); % Supply same shocks regardless of parameters
        catch ME
            disp('Error encountered in likelihood computation. Message:');
            disp(ME.message);
        end
        
    end

    % Print progress
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_lik, '/', lik_numgrid);

end

%% Save results

mkdir(save_folder);
save_mat(fullfile(save_folder, model_name));

delete(poolobj);


%% Plot

Z = reshape(lik_all(:,1),param_vals_num,[]);
Z = Z - max(Z(:));
% contour(param_vals(:,2),param_vals(:,1),Z,0:-10:-200);
contour(param_vals(:,2),param_vals(:,1),Z,0:-1:-20);
[~,max_ind]=max(lik_all(:,1));
hold on;
plot(params_truth(2),params_truth(1),'xr','MarkerSize',10);
plot(lik_grid(max_ind,2),lik_grid(max_ind,1),'ok','MarkerSize',10);
hold off;
xlabel(param_names{2});
ylabel(param_names{1});

% xlim(params_truth(2)*[0.9 1.1]);
% ylim(params_truth(1)*[0.9 1.1]);