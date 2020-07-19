clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim', 'auxiliary_functions/mcmc');


%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation
likelihood_type = 1; % =1: Macro + full-info micro; =2: macro + full-info micro, no truncation; =3: macro + moments micro


% Model/data settings
T = 100;                                % Number of periods of simulated macro data
ts_micro = 10:10:T;                        % Time periods where we observe micro data
N_micro = 1e2;                             % Number of households per non-missing time period

% Parameter transformation
transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2)) -exp(x(3))]; % Function mapping transformed parameters into parameters of interest
param_to_transf = @(x) [log(x(1)/(1-x(1))) log(x(2)) log(-x(3))];  % Function mapping parameters of interest into transformed parameters

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters

% Initialization settings
init_type = 3; % 1: start from true parameter values
               % 2: start from other values
               % 3: find posterior mode via optimization
if init_type == 3
%     optim_grid = combvec(linspace(0.8,0.99,5),linspace(0.01,0.05,5),linspace(-.5,-0.1,5))';    % Optimization grid
    [aux1, aux2, aux3] = meshgrid(linspace(0.8,0.99,5),linspace(0.01,0.05,5),linspace(-.5,-0.1,5)); % same purpose, no need of Deep Learning Toolbox
    optim_grid = [aux1(:), aux2(:), aux3(:)];
end

% MCMC settings
mcmc_num_draws = 9000;                  % Number of MCMC steps (total)
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [50 100 200];          % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = false;                 % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                  % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)
mcmc_filename = ['mcmc_N' num2str(N_micro) '_liktype' num2str(likelihood_type)...
    '_inittype' num2str(init_type) '.mat'];              % File name of MCMC output

% for adaptive RWMH
mcmc_c = 0.55;
mcmc_ar_tg = 0.3;
p_adapt = .95;

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 201807271;                    % Random number generator seed for initial simulation
rng('default');

%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP ttau mu_l ssigmaMeas;
	
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

% Distribution of indv params log(lambda_i) ~ N(mu_l,-2*mu_l), 
% so lambda_i > 0 and E(lambda_i) = 1
mu_l = -.25; % Roughly calibrated to Piketty, Saez & Zucman (QJE 2018) Table I, log of P90-P20 ratio, post-tax

% Measurement error std. dev. of observed log output
ssigmaMeas = 0.02;


%% Set approximation parameters

global nEpsilon nAssets nAssetsFine nAssetsQuadrature ...
	nMeasure maxIterations tolerance dampening splineOpt displayOpt;

% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% MUST BE SET TO 0

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

cd('./auxiliary_functions/dynare');
delete steady_vars.mat;
saveParameters;


%% Initial Dynare run

rng(rng_seed);
dynare firstOrderDynamics_polynomials noclearall nopathchange; % Run Dynare once to process model file


%% Simulate data

if is_data_gen == 1

    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    for i_Epsilon = 1:nEpsilon
        for i_Measure = 1:nMeasure
            sim_struct.(sprintf('%s%d%d', 'smpl_m', i_Epsilon, i_Measure)) = nan(T,1); 
                % Set sample moments to missing everywhere
        end
    end
    % Add measurement error to aggregate output
    sim_struct.logAggregateOutput_noerror = sim_struct.logAggregateOutput;
    sim_struct.logAggregateOutput = sim_struct.logAggregateOutput + ssigmaMeas*randn(T,1);
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data

    % draw normalized individual incomes
    simul_data_micro_aux = simulate_micro_aux(sim_struct, ts_micro, N_micro);

    % draw individual productivities and incomes
    simul_data_micro = simulate_micro(simul_data_micro_aux);
    save('simul_data_micro.mat','simul_data_micro');

    % Compute cross-sectional moments from micro data
    sim_struct_moments = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro);
    save('simul_moments.mat', '-struct', 'sim_struct_moments');
    
    micro_mom_order = 2;
    for im=micro_mom_order+1:nMeasure
        % Set higher-order sample moments to missing
        sim_struct_moments.(sprintf('%s%d', 'smpl_m1', im)) = nan(T,1);
        sim_struct_moments.(sprintf('%s%d', 'smpl_m2', im)) = nan(T,1);
    end
    save('simul_moments2.mat', '-struct', 'sim_struct_moments');

else

    % Load previous data
%         load('simul.mat')
    load('simul_data_micro.mat');

end

%% Optimize over grid to find approximate mode

if init_type == 3
    
    optim_numgrid = size(optim_grid,1);
    optim_logpost = nan(optim_numgrid,1);
    
    disp('Optimization...');
    
    for i_optim=1:optim_numgrid % Cycle through grid points
        
        the_param = optim_grid(i_optim,:);
        fprintf(['%s' repmat('%6.4f ',1,length(the_param)),'%s\n'], 'current  [bbeta,ssigmaMeas,mu_l] = [',...
        the_param,']');
        the_aux = num2cell(the_param);
        [bbeta,ssigmaMeas,mu_l] = deal(the_aux{:});
        
        try
            saveParameters;         % Save parameter values to files
            setDynareParameters;    % Update Dynare parameters in model struct
            compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
            compute_meas_err;       % Update measurement error var-cov matrix
        
            the_loglike = aux_compute_likel(simul_data_micro, ts_micro, ...
                                   num_smooth_draws, num_interp, num_burnin_periods, ...
                                   M_, oo_, options_, ...
                                   likelihood_type);
            the_logprior = prior_logdens_transf(param_to_transf(the_param));
            optim_logpost(i_optim) = the_loglike + the_logprior;
        catch ME
            disp('Error encountered. Message:');
            disp(ME.message);
        end
    
        % Print progress
        fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_optim, '/', optim_numgrid);
        
    end
    
    % Find maximum of log posterior on grid
    [~,optim_maxind] = max(optim_logpost);
    optim_maxparam = optim_grid(optim_maxind,:);
    
    fprintf(['%s' repmat('%6.4f ',1,length(optim_maxparam)),'%s\n'], 'approximate mode  [bbeta,ssigmaMeas,mu_l] = [',...
        optim_maxparam,']');

end


%% MCMC

if init_type == 1
    init_val = [bbeta,ssigmaMeas,mu_l];
elseif init_type == 2
    init_val = [.9 .06 -1];
elseif init_type == 3
    init_val = optim_maxparam;
end

curr_draw = param_to_transf(init_val);        % Initial draw
post_draws = nan(mcmc_num_draws,length(curr_draw));
accepts = zeros(mcmc_num_draws,1);

curr_logpost = -Inf;
loglikes_prop = nan(mcmc_num_draws,1);
loglikes_prop_macro = nan(mcmc_num_draws,1);
loglikes_prop_micro = nan(mcmc_num_draws,1);

the_stepsize = mcmc_stepsize_init;      % Initial RWMH step size
the_stepsize_iter = 1;
the_chol = eye(length(curr_draw));      % Initial RWMH proposal var-cov matrix

disp('MCMC...');
timer_mcmc = tic;

if likelihood_type == 1
    delete(gcp('nocreate'));
    
    if contains(pwd,'u050')
        parpool('local',12);
    else
        parpool;
    end
end

for i_mcmc=1:mcmc_num_draws % For each MCMC step...

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'current  [bbeta,ssigmaMeas,mu_l] = [',...
        transf_to_param(curr_draw),']');
    
    % Proposed draw (modified to always start with initial draw)
    [prop_draw,is_adapt] = rwmh_propose(curr_draw, (i_mcmc>1)*the_stepsize, the_chol, p_adapt, mcmc_stepsize_init); % Proposal
    
    % Set new parameters
    the_transf = num2cell(transf_to_param(prop_draw));
    [bbeta,ssigmaMeas,mu_l] = deal(the_transf{:});

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'proposed [bbeta,ssigmaMeas,mu_l] = [',...
        [bbeta,ssigmaMeas,mu_l],']');
    
    try

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
        compute_meas_err;       % Update measurement error var-cov matrix

        [loglikes_prop(i_mcmc), loglikes_prop_macro(i_mcmc), loglikes_prop_micro(i_mcmc)] = ...
            aux_compute_likel(simul_data_micro, ts_micro, ...
            num_smooth_draws, num_interp, num_burnin_periods, ...
            M_, oo_, options_, ...
            likelihood_type);

        % Log prior density of proposal
        logprior_prop = prior_logdens_transf(prop_draw);

        % Accept/reject
        [curr_draw, curr_logpost, accepts(i_mcmc), the_log_ar] = rwmh_accrej(curr_draw, prop_draw, curr_logpost, logprior_prop+loglikes_prop(i_mcmc));

        % Adapt proposal step size
        if is_adapt
            [the_stepsize, the_stepsize_iter] = adapt_stepsize(the_stepsize, the_stepsize_iter, i_mcmc, the_log_ar, mcmc_c, mcmc_ar_tg);
        end
        
    catch ME
        
        disp('Error encountered. Message:');
        disp(ME.message);
        
    end
    
    % Store
    post_draws(i_mcmc,:) = transf_to_param(curr_draw);
    
    % Print acceptance rate
    fprintf('%s%5.1f%s\n', 'Accept. rate last 100: ', 100*mean(accepts(max(i_mcmc-99,1):i_mcmc)), '%');
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_mcmc, '/', mcmc_num_draws);
    
    % Adapt proposal covariance matrix
    [the_chol, the_stepsize_iter] = adapt_cov(the_chol, the_stepsize_iter, mcmc_adapt_iter, i_mcmc, post_draws, mcmc_adapt_diag, mcmc_adapt_param);
    
    % Save middle steps in case reach time limit on the server
    if mod(i_mcmc,1000) == 0
        save(mcmc_filename);
    end
    
end

if likelihood_type == 1
    delete(gcp('nocreate'));
end

mcmc_elapsed = toc(timer_mcmc);
fprintf('%s%8.2f\n', 'MCMC done. Elapsed minutes: ', mcmc_elapsed/60);

cd('../../');

save(mcmc_filename);
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');

%% Auxiliary function

function [the_loglike, the_loglike_macro, the_loglike_micro] = ...
         aux_compute_likel(simul_data_micro, ts_micro, ...
                           num_smooth_draws, num_interp, num_burnin_periods, ...
                           M_, oo_, options_, ...
                           likelihood_type)

        % Log likelihood of proposal
        switch likelihood_type
            case 1 % Macro + full info micro
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                    num_smooth_draws, num_interp, num_burnin_periods, ...
                    M_, oo_, options_);
            case 2 % Macro + moments w/ SS meas. err.
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute('simul_moments.mat', [], ts_micro, ...
                    num_smooth_draws, -Inf, num_burnin_periods, ...
                    M_, oo_, options_);
            case 3  % Macro + moments w/ SS meas. err. (observe up to 2nd moment)
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute('simul_moments2.mat', [], ts_micro, ...
                    0, num_interp, num_burnin_periods, ...
                    M_, oo_, options_);
        end
        
end