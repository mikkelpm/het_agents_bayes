clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');


%% Settings

% Decide what to do
is_data_gen = 2; % whether simulate data:  
                 % 0: no simulation; 
                 % 1: simulation (with indv params, based on simulation without indv params)
                 % 2: simulation (with indv params, start from scratch)
is_profile = 0; %whether run profiler for execution time

% Model/data settings
T = 200;                                % Number of periods of simulated macro data
ts_hh = 20:20:T;                      % Time periods where we observe micro data
N_hh = 1e3;                              % Number of households per non-missing time period

% Parameter transformation
transf_to_param = @(x) [exp(x(1)) exp(x(2)) -exp(x(3))]; % Function mapping transformed parameters into parameters of interest

% Prior
prior_logdens_transf = @(x) sum(x) + ((-1-1)*x(1)-4/exp(x(1)));    % Log prior density of transformed parameters
% prior_logdens_log = @(x) sum(x)-2*log(1+exp(x(3)));    % Log prior density of log(beta)
prior_init_transf = @() [log(1) log(1) log(.5)];  % Distribution of initial transformed parameter draw
% prior_init_log = @() [log(1) log(1) log(.895)-log(1-.895) log(.014) log(.5)];  % Distribution of initial log(beta) draw

% MCMC settings
mcmc_num_draws = 500;                   % Number of MCMC steps (total)
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [20 50 100];         % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
mcmc_adapt_diag = true;                    % =true: Adapt only to posterior std devs of parameters, =false: adapt to full var/cov matrix
mcmc_adapt_param = 10;                    % Shrinkage parameter for adapting to var/cov matrix (higher values: more shrinkage)
% mcmc_num_adapt = 50;                    % Number of initial steps with larger (but gradually decreasing) step size
% mcmc_stepsizes = [1e-3 1e-2];           % MCMC step size after initial phase
% mcmc_stepsizes_init = 10*mcmc_stepsizes; % Initial MCMC steps size (gradually decreased during initial phase)
mcmc_smooth_draws = 500;                % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
mcmc_filename = 'mcmc.mat';      % File name of MCMC output

% for adaptive RWMH
mcmc_c = 0.55;
mcmc_ar_tg = 0.3;

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20180709;                    % Random number generator seed for initial simulation

% Profiler save settings
tag_date = datestr(now,'yyyymmdd');


%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP ttau mu_l;
	
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

% Distribution of indv params log(lambda_i) ~ N(-1/2,1), 
% so lambda_i > 0 and E(lambda_i) = 1
mu_l = -.5;


%% Set approximation parameters

global nEpsilon nAssets nAssetsFine nAssetsQuadrature ...
	nMeasure maxIterations tolerance dampening splineOpt displayOpt;

% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% if splineOpt = 1, use splines to approximate savings policy; if splineOpt = 0, use polynomials
				% to approximate conditional expectation function

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
dynare firstOrderDynamics_polynomials noclearall;                   % Run Dynare once to process model file


%% Simulate data

if is_data_gen == 0
    
    load('simul.mat')
    load('simul_data_hh_indv_param.mat');
    
elseif is_data_gen == 1
    
    load('simul.mat')
    load('simul_data_hh.mat');
    
    % draw individual productivities and incomes
    simul_data_hh_indv_param = simulate_hh_indv_param(simul_data_hh);
    save('simul_data_hh_indv_param.mat','simul_data_hh_indv_param');
    
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

curr_loglike = -Inf;
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

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'current  [ssigma,uDuration,mu_l] = [',...
        transf_to_param(curr_draw),']');
%     fprintf(['%s' repmat('%6.4f',1,5),'%s\n'], 'current  [ssigma,uDuration,rrhoTFP,ssigmaTFP,mu_l] = [',...
%         [exp(curr_log([1 2])) 1/(1+exp(-curr_log(3))) exp(curr_log(4)) -exp(curr_log(5))],']');
    
    % Proposed draw (modified to always start with initial draw)
    prop_draw = curr_draw + (i_mcmc>1)*the_stepsize*randn(size(curr_draw))*the_chol; % Proposal
    
    % Set new parameters
    the_transf = num2cell(transf_to_param(prop_draw));
    [ssigma,uDuration,mu_l] = deal(the_transf{:});

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'proposed [ssigma,uDuration,mu_l] = [',...
        [ssigma,uDuration,mu_l],']');
%     fprintf(['%s' repmat('%6.4f',1,5),'%s\n'], 'poposed  [ssigma,uDuration,rrhoTFP,ssigmaTFP,mu_l] = [',...
%         [ssigma,uDuration,rrhoTFP,ssigmaTFP,mu_l],']');
    
    try
        
        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state once and for all
        
        % Log likelihood of proposal
        [loglikes_prop(i_mcmc), loglikes_prop_macro(i_mcmc), loglikes_prop_hh(i_mcmc)] = ...
            loglike_compute_indv_param('simul.mat', simul_data_hh_indv_param, ts_hh, mcmc_smooth_draws, num_burnin_periods, M_, oo_, options_);

        % Log prior ratio
        logratio_prior = prior_logdens_transf(prop_draw)-prior_logdens_transf(curr_draw);

        log_ar = min(loglikes_prop(i_mcmc)-curr_loglike+logratio_prior,0);
        % Accept/reject
        if log_ar > log(rand())
            fprintf('%s\n', 'Accepted!');
            curr_draw = prop_draw;
            curr_loglike = loglikes_prop(i_mcmc);
            accepts(i_mcmc) = 1;
        end

        % Adaptive step size, cf. Atchade and Rosenthal (2005), Section 4.1
        the_stepsize = exp(log(the_stepsize)+(i_mcmc>1)*the_stepsize_iter^(-mcmc_c)*(exp(log_ar)-mcmc_ar_tg)); 
        the_stepsize = min(max(the_stepsize,1e-10),1e10);
        the_stepsize_iter = the_stepsize_iter + 1;
        fprintf('%s%8.6f\n', 'New step size: ', the_stepsize);
    
    catch ME
        
        disp('Error encountered. Message:');
        disp(ME.message);
        
    end
    
    % Store
    post_draws(i_mcmc,:) = transf_to_param(curr_draw);
    
    % Print acceptance rate
    fprintf('%s%5.1f%s\n', 'Accept. rate last 100: ', 100*mean(accepts(max(i_mcmc-99,1):i_mcmc)), '%');
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_mcmc, '/', mcmc_num_draws);
    
    % Update RWMH proposal var-cov matrix
    the_indx = find(mcmc_adapt_iter==i_mcmc,1); % Find current iteration index in list of adaptation iterations
    if ~isempty(the_indx) && the_indx>1 % If in list, but not first...
        the_start = 1 + (the_indx>1)*mcmc_adapt_iter(the_indx-1); % Start of current adaptation window
        the_cov = cov(post_draws(the_start:i_mcmc,:)); % Var-cov matrix over current adaptation window
        the_cov = the_cov/(mean(sqrt(diag(the_cov)))^2); % Normalize average std deviation
        the_n = i_mcmc-the_start+1;
        the_cov_shrink = (the_n*the_cov + mcmc_adapt_param*eye(length(curr_draw))) ...
                         /(the_n+mcmc_adapt_param); % Shrink var-cov matrix toward identity
        if mcmc_adapt_diag % If adapt only to std devs...
            the_chol = diag(sqrt(diag(the_cov_shrink)));
        else % If adapt full var-cov matrix...
            the_chol = chol(the_cov_shrink);
        end
        disp('New RWMH proposal var-cov matrix:');
        disp(the_chol'*the_chol);
        disp('Square root of diagonal');
        disp(sqrt(diag(the_chol'*the_chol)));
        the_stepsize_iter = 1; % Reset stepsize adaptation
    end
    
end

delete(poolobj);

mcmc_elapsed = toc(timer_mcmc);
fprintf('%s%8.2f\n', 'MCMC done. Elapsed minutes: ', mcmc_elapsed/60);

cd('../../');

save(mcmc_filename); %, 'mcmc_*', 'post_draws', 'accepts', 'loglikes_prop*');

if is_profile
    profsave(profile('info'),['profile_results_' tag_date]);
end