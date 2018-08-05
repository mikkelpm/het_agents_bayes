clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim', 'auxiliary_functions/mcmc');


%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation

% Model/data settings
T = 100;                                % Number of periods of simulated macro data
ts_hh = 10:10:T;                        % Time periods where we observe micro data
N_hh = 1e3;                             % Number of households per non-missing time period

% Parameter transformation
transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2)) -exp(x(3))]; % Function mapping transformed parameters into parameters of interest

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters
prior_init_transf = @() [log(0.96)-log(1-0.96) log(0.02) log(.25)];  % Distribution of initial log(beta) draw

% MCMC settings
mcmc_num_draws = 1000;                  % Number of MCMC steps (total)
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [50 100 200];          % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = false;                 % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                  % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)
mcmc_filename = 'mcmc.mat';             % File name of MCMC output

% for adaptive RWMH
mcmc_c = 0.55;
mcmc_ar_tg = 0.3;

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 201807271;                    % Random number generator seed for initial simulation

% Profiler save settings
tag_date = datestr(now,'yyyymmdd');


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

if is_data_gen == 0
    
    % Load previous data
    load('simul.mat')
    load('simul_data_hh_indv_param.mat');
    
else
    
    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw normalized individual incomes
    simul_data_hh = simulate_hh(sim_struct, ts_hh, N_hh);
    save('simul_data_hh.mat','simul_data_hh');
    
    % draw individual productivities and incomes
    simul_data_hh_indv_param = simulate_hh_indv_param(simul_data_hh);
    save('simul_data_hh_indv_param.mat','simul_data_hh_indv_param');
    
end


%% MCMC

curr_draw = prior_init_transf();        % Initial draw
post_draws = nan(mcmc_num_draws,length(curr_draw));
accepts = zeros(mcmc_num_draws,1);

curr_logpost = -Inf;
loglikes_prop = nan(mcmc_num_draws,1);
loglikes_prop_macro = nan(mcmc_num_draws,1);
loglikes_prop_hh = nan(mcmc_num_draws,1);

the_stepsize = mcmc_stepsize_init;      % Initial RWMH step size
the_stepsize_iter = 1;
the_chol = eye(length(curr_draw));      % Initial RWMH proposal var-cov matrix

disp('MCMC...');
timer_mcmc = tic;

poolobj = parpool;

for i_mcmc=1:mcmc_num_draws % For each MCMC step...

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'current  [bbeta,ssigmaMeas,mu_l] = [',...
        transf_to_param(curr_draw),']');
    
    % Proposed draw (modified to always start with initial draw)
    prop_draw = rwmh_propose(curr_draw, (i_mcmc>1)*the_stepsize, the_chol); % Proposal
    
    % Set new parameters
    the_transf = num2cell(transf_to_param(prop_draw));
    [bbeta,ssigmaMeas,mu_l] = deal(the_transf{:});

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'proposed [bbeta,ssigmaMeas,mu_l] = [',...
        [bbeta,ssigmaMeas,mu_l],']');
    
    try
        
        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state once and for all
        
        % Log likelihood of proposal
        [loglikes_prop(i_mcmc), loglikes_prop_macro(i_mcmc), loglikes_prop_hh(i_mcmc)] = ...
            loglike_compute_indv_param('simul.mat', simul_data_hh_indv_param, ts_hh, ...
                                       num_smooth_draws, num_interp, num_burnin_periods, ...
                                       M_, oo_, options_);

        % Log prior density of proposal
        logprior_prop = prior_logdens_transf(prop_draw);

        % Accept/reject
        [curr_draw, curr_logpost, accepts(i_mcmc), the_log_ar] = rwmh_accrej(curr_draw, prop_draw, curr_logpost, logprior_prop+loglikes_prop(i_mcmc));

        % Adapt proposal step size
        [the_stepsize, the_stepsize_iter] = adapt_stepsize(the_stepsize, the_stepsize_iter, i_mcmc, the_log_ar, mcmc_c, mcmc_ar_tg);
    
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
    
end

delete(poolobj);

mcmc_elapsed = toc(timer_mcmc);
fprintf('%s%8.2f\n', 'MCMC done. Elapsed minutes: ', mcmc_elapsed/60);

cd('../../');

save(mcmc_filename);
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');