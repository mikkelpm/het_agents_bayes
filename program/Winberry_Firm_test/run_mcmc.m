clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/mcmc', 'auxiliary_functions/sim');


%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation

% Model/data settings
T = 50;                                % Number of periods of simulated macro data
ts_micro = 10:10:T;                        % Time periods where we observe micro data
N_micro = 1e3;                             % Number of micro entities per non-missing time period

% Parameter transformation
transf_to_param = @(x) [1/(1+exp(-x(1))) exp(x(2:end))]; % Function mapping transformed parameters into parameters of interest

% Prior
prior_logdens_transf = @(x) sum(x) - 2*log(1+exp(x(1)));    % Log prior density of transformed parameters
prior_init_transf = @() [log(0.53)-log(1-0.53) log(0.0364) log(.011) log(.0083)];  % Distribution of initial transformed draw

% MCMC settings
mcmc_num_iter = 1e4;                  % Number of MCMC steps (total)
mcmc_thin = 2;                         % Store every X draws
mcmc_stepsize_init = 1e-2;              % Initial MCMC step size
mcmc_adapt_iter = [5 200 500 1000];          % Iterations at which to update the variance/covariance matrix for RWMH proposal; first iteration in list is start of adaptation phase
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
rng_seed = 201809251;                    % Random number generator seed for initial simulation

%% Set economic parameters 

global ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ  cchi

% Technology
ttheta 			= .256;								% capital coefficient
nnu 				= .64;								% labor coefficient
ddelta 			= .085;								% depreciation (annual)
rrhoProd 		= .53; %.859; 								% persistence of idiosyncratic shocks (annual)
ssigmaProd 	= .0364; %.022;								% SD innovations of idiosycnratic shocks (annual)
aaUpper 		= .011; 								% no fixed cost region upper bound
aaLower 		= -.011;								% no fixed cost region lower bound
ppsiCapital 	= .0083;							% upper bound on fixed adjustment cost draws
% ppsiCapital     = 1e-5; 

% Preferences
bbeta 			= .961;								% discount factor (annual)
ssigma 			= 1;									% coefficient of relative risk aversion
pphi 				= 1 / 1e5;							% inverse Frisch elasticity of labor supply
nSS 				= 1 / 3;								% hours worked in steady state
cchi				= 1;									% labor disutility (will be calibrated in Dynare's steady state file to ensure labor supply = nSS)

% Aggregate shocks
rrhoTFP			= 0.859;							% persistence of aggregate TFP (annual)
ssigmaTFP		= .014;								% SD of innovations of aggregate TFP (annual)
rrhoQ				= 0.859;							% persistence of aggregate investment-specific shock (annual)
ssigmaQ		= .014;								% SD of innovations of aggregate investment-specific shock (annual)
corrTFPQ		= 0;									% loading on TFP shock in evolution of investment-specific shock

%% Set approximation parameters

global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
	maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
	nCapitalQuadrature kRepSS wRepSS

% Order of approximation of value function
nProd 			= 3;										% order of polynomials in productivity
nCapital 		= 5;										% order of polynomials in capital
nState 			= nProd * nCapital;					% total number of coefficients

% Shocks 
nShocks 		= 3;										% order of Gauss-Hermite quadrature over idiosyncratic shocks

% Finer grid for analyzing policy functions and computing histogram
nProdFine 		= 60;
nCapitalFine 	= 40;
nStateFine 	= nProdFine * nCapitalFine;

% Iteration on value function
maxIterations	= 100;
tolerance 		= 1e-6;
acc 				= 500;									% number of iterations in "Howard improvement step"
dampening 	= 0;										% weight on old iteration in updating step

% Approximation of distribution
nMeasure 				= 2;							% order of polynomial approximating distribution
nProdQuadrature 		= 8; 							% number of quadrature points in productivity dimension
nCapitalQuadrature 	= 10;						% number of quadrature points in capital dimension
nStateQuadrature 		= nProdQuadrature * nCapitalQuadrature;
nMeasureCoefficients 	= (nMeasure * (nMeasure + 1)) / 2 + nMeasure;

%% Save parameters

cd('./auxiliary_functions/dynare');
delete steady_vars.mat;
saveParameters;


%% Initial Dynare run

rng(rng_seed);
dynare dynamicModel noclearall nopathchange; % Run Dynare once to process model file


%% Simulate data

if is_data_gen == 0
    
    % Load previous data
    load('simul.mat')
    load('simul_data_micro.mat');
%     load('simul_data_micro_indv_param.mat');
    
else
    
    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw micro data
    simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, num_interp);
    save('simul_data_micro.mat','simul_data_micro');
    
%     % draw individual productivities and incomes
%     simul_data_micro_indv_param = simulate_micro_indv_param(simul_data_micro);
%     save('simul_data_micro_indv_param.mat','simul_data_micro_indv_param');
    
end


%% MCMC

curr_draw = prior_init_transf();        % Initial draw

accepts = zeros(mcmc_num_iter,1);

mcmc_num_draws = floor(mcmc_num_iter/mcmc_thin); % Number of stored draws
post_draws = nan(mcmc_num_draws,length(curr_draw));

curr_logpost = -Inf;
loglikes_prop = nan(mcmc_num_draws,1);
loglikes_prop_macro = nan(mcmc_num_draws,1);
loglikes_prop_micro = nan(mcmc_num_draws,1);

the_stepsize = mcmc_stepsize_init;      % Initial RWMH step size
the_stepsize_iter = 1;
the_chol = eye(length(curr_draw));      % Initial RWMH proposal var-cov matrix

disp('MCMC...');
timer_mcmc = tic;

poolobj = parpool;

for i_mcmc=1:mcmc_num_iter % For each MCMC step...
    
    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'current  [rrhoProd,ssigmaProd,aaUpper,ppsiCapital] = [',...
        transf_to_param(curr_draw),']');
    
    % Proposed draw (modified to always start with initial draw)
    prop_draw = rwmh_propose(curr_draw, (i_mcmc>1)*the_stepsize, the_chol); % Proposal
    
    % Set new parameters
    the_transf = num2cell(transf_to_param(prop_draw));
    [rrhoProd,ssigmaProd,aaUpper,ppsiCapital] = deal(the_transf{:});
    aaLower = -aaUpper;

    fprintf(['%s' repmat('%6.4f ',1,length(curr_draw)),'%s\n'], 'proposed [rrhoProd,ssigmaProd,aaUpper,ppsiCapital] = [',...
        [rrhoProd,ssigmaProd,aaUpper,ppsiCapital],']');
    
    the_loglike_prop = [];
    the_loglike_prop_macro = [];
    the_loglike_prop_micro = [];
    
    try
        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics

        % Log likelihood of proposal
        [the_loglike_prop, the_loglike_prop_macro, the_loglike_prop_micro] = ...
            loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                                       num_smooth_draws, num_interp, num_burnin_periods, ...
                                       M_, oo_, options_);
        
        % Log prior density of proposal
        logprior_prop = prior_logdens_transf(prop_draw);

        % Accept/reject
        [curr_draw, curr_logpost, accepts(i_mcmc), the_log_ar] = rwmh_accrej(curr_draw, prop_draw, curr_logpost, logprior_prop+the_loglike_prop);

        % Adapt proposal step size
        [the_stepsize, the_stepsize_iter] = adapt_stepsize(the_stepsize, the_stepsize_iter, i_mcmc, the_log_ar, mcmc_c, mcmc_ar_tg);
    
    catch ME
        
        disp('Error encountered. Message:');
        disp(ME.message);
        
    end
    
    % Store (with thinning)
    if mod(i_mcmc, mcmc_thin)==0
        post_draws(i_mcmc/mcmc_thin,:) = transf_to_param(curr_draw);
        loglikes_prop(i_mcmc/mcmc_thin) = the_loglike_prop;
        loglikes_prop_macro(i_mcmc/mcmc_thin) = the_loglike_prop_macro;
        loglikes_prop_micro(i_mcmc/mcmc_thin) = the_loglike_prop_micro;
    end
    
    % Print acceptance rate
    fprintf('%s%5.1f%s\n', 'Accept. rate last 100: ', 100*mean(accepts(max(i_mcmc-99,1):i_mcmc)), '%');
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_mcmc, '/', mcmc_num_iter);
    
    % Adapt proposal covariance matrix
    [the_chol, the_stepsize_iter] = adapt_cov(the_chol, the_stepsize_iter, mcmc_adapt_iter, i_mcmc, post_draws, mcmc_thin, mcmc_adapt_diag, mcmc_adapt_param);
    
end

delete(poolobj);

mcmc_elapsed = toc(timer_mcmc);
fprintf('%s%8.2f\n', 'MCMC done. Elapsed minutes: ', mcmc_elapsed/60);

cd('../../');

save(mcmc_filename);
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');