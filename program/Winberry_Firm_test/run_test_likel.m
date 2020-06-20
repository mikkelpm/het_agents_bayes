clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');


%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation
is_profile = 0; %whether run profiler for execution time
if is_profile == 1
    profile on
end

% Model/data settings
T = 50;                                % Number of periods of simulated macro data
ts_micro = 10:10:T;                        % Time periods where we observe micro data
N_micro = 1e2;                             % Number of micro entities per non-missing time period

% Parameter values to check (prod dynamics)
param_plot = 'rrhoProd'; %'ssigmaProd'; %'ppsiCapital'; %'aaUpper'; % Parameter to vary when plotting likelihood
param_vals_mult = unique([1 linspace(0.75,1.25,50)]); % Multiples of true parameter to compute

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 202006171;                   % Random number generator seed for initial simulation


%% Set economic parameters 

global ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ  cchi

% Technology
ttheta 			= .256;								% capital coefficient
nnu 				= .64;								% labor coefficient
ddelta 			= .085;								% depreciation (annual)
rrhoProd 		= .53; %.859; 								% persistence of idiosyncratic shocks (annual)
ssigmaProd 	= .0364;  %.022;								% SD innovations of idiosycnratic shocks (annual)
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
% return;


%% Simulate data

if is_data_gen == 0
    
    % Load previous data
%     load('simul.mat')
    load('simul_data_micro.mat');
%     load('simul_data_micro_indv_param.mat');
    
else
    
    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    for j=1:nMeasureCoefficients
        sim_struct.(sprintf('%s%d', 'smpl_m', j)) = nan(T,1); % Set sample moments to missing everywhere
    end
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw micro data
    simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, num_interp);
    save('simul_data_micro.mat','simul_data_micro');
    
    % Compute cross-sectional moments from micro data
    sim_struct_moments = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro);
    save('simul_moments.mat', '-struct', 'sim_struct_moments');
    
%     % draw individual productivities and incomes
%     simul_data_micro_indv_param = simulate_micro_indv_param(simul_data_micro);
%     save('simul_data_micro_indv_param.mat','simul_data_micro_indv_param');
    
end


%% Compute likelihood

eval(['param_true = ' param_plot ';']);
param_vals = param_true*param_vals_mult;
loglikes = nan(length(param_vals),3);
loglikes_macro = nan(length(param_vals),3);
loglikes_micro = nan(length(param_vals),3);

disp('Computing likelihood...');
timer_likelihood = tic;

if contains(pwd,'Laura')
    poolobj = parpool(2);
else
    poolobj = parpool;
end

for iter_i=1:length(param_vals) % For each parameter value...
        
    % Set new parameter
    eval([param_plot ' = param_vals(iter_i);']);
    fprintf('%d%s%d%s%s%s%6.4f\n', iter_i, '/', length(param_vals), ': ', param_plot, ' = ', eval(param_plot));

    saveParameters;         % Save parameter values to files
    setDynareParameters;    % Update Dynare parameters in model struct
    compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
    compute_meas_err;       % Update measurement error var-cov matrix

    % Log likelihood

    % Macro + full info micro
    [loglikes(iter_i,1), loglikes_macro(iter_i,1), loglikes_micro(iter_i,1)] = ...
        loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                                   num_smooth_draws, num_interp, num_burnin_periods, ...
                                   M_, oo_, options_);
    % Macro + moments w/ SS meas. err.
    [loglikes(iter_i,2), loglikes_macro(iter_i,2), loglikes_micro(iter_i,2)] = ...
        loglike_compute('simul_moments.mat', [], ts_micro, ...
                                   0, num_interp, num_burnin_periods, ...
                                   M_, oo_, options_);
    % Macro + moments w/o meas. err.
    M_.H(3:end,3:end) = 1e-8*eye(5); % Only a little bit of meas. err. to avoid singularity
    [loglikes(iter_i,3), loglikes_macro(iter_i,3), loglikes_micro(iter_i,3)] = ...
        loglike_compute('simul_moments.mat', [], ts_micro, ...
                                   0, num_interp, num_burnin_periods, ...
                                   M_, oo_, options_);
    
end

delete(poolobj);

likelihood_elapsed = toc(timer_likelihood);
fprintf('%s%8.2f\n', 'Done. Elapsed minutes: ', likelihood_elapsed/60);

cd('../../');
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');


%% Plot

loglikes_plot = [loglikes(:,1:2) loglikes_macro(:,1)];
plot(param_vals,loglikes_plot-max(loglikes_plot,[],1));
the_ylim = ylim;
line(param_true*ones(1,2), the_ylim, 'Color', 'k', 'LineStyle', ':');
ylim(the_ylim);
legend({'FI micro', 'moments micro', 'macro'}, 'Location', 'SouthWest');


if is_profile
    profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);
end
