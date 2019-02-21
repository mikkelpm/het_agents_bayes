% Shell which declares parameters and calls Dynare to solve model using
% first order approximation of aggregate dynamics
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');

% profile on

cd('./auxiliary_functions/dynare');


%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation

% Model/data settings
T = 200;                        % Number of periods of simulated macro data
ts_micro = 20:20:T;             % Time periods where we observe micro data
N_micro = 1e3;                  % Number of households per non-missing time period

% Parameter values to check
param1_vals = [0.98 0.99 0.999];
param2_vals = [50 100 200];

% Likelihood settings
num_smooth_draws = 500;         % Number of draws from the smoothing distribution (for unbiased likelihood estimate)

% Numerical settings
num_burnin_periods = 100;       % Number of burn-in periods for simulations
rng_seed = 201902211;            % Random number generator seed for initial simulation


%% Define parameters

% ECONOMIC PARAMETERS

% Preferences
bbeta = .99; % Discount factor (quarterly)
ppsi = 1; % Coefficient on labor disutility
nnu = 1; % Inverse Frisch elasticity, MUST be 1 at the moment!

% Borrowing limit
bbBar = -1;

% Production
eepsilon = 5;

% Taxes
ttau = 0.3;
vvarthetaB = -0.233;
vvarthetaT = 0.06;

% Idiosyncratic productivity
vzGrid = [0.5;1;1.5];
mzTransition = [0.8 0.2 0;
                0.1 0.8 0.1;
                0   0.2 0.8]; % (i,j) element: P(z'=z_j | z=z_i)

% Aggregate productivity
A_SS = 3; % Steady state aggr productivity level

% Equity shares
% vShareGrid = [0; 1; 2]; % Profit shares for each household type
% vShareFraction = [1/3; 1/3; 1/3]; % Fractions of each household type
vShareGrid = 1;
vShareFraction = 1;

% Aggregate shocks
rrhoTFP = .95; % Autocorrelation of TFP
ssigmaTFP = .005; % Std dev of innovation to TFP shock
rrhoMP = 0.9; % Autocorrelation of monetary policy shock
ssigmaMP = 0.0025; % Std dev of innovation to monetary policy shock

% Dynamic parameters (KMV (2018))
ttheta = 100; % Rotemberg price stickiness parameter
pphi = 1.25; % Taylor rule coefficient on inflation

% APPROXIMATION PARAMETERS

% Whether to print out results from steady state computation
displayOpt = 'off';       % 'iter-detailed' or 'off'

% Order of approximation
nAssets = 25; % number of polynomials in polynomial approximation

% Finer grid for analyzing policy functions
nAssetsFine = 100;

% Approximation of distribution
nMeasure = 3; % Polynomial order of parametric density approximation
nAssetsQuadrature = 8; % Number of quadrature points for assets

% Steady state
tolerance_SS_root = 0.001; % Numerical tolerance for root finding
tolerance_SS_invhist = 1e-12; % Numerical tolerance for invariant distribution of histogram approach
maxIterations = 2e4; % Max no. of iterations in inner steady state loops
tolerance = 1e-5; % Numerical tolerance for inner steady state loops
dampening = .5; % Dampening factor in steady state iterations (higher value -> more dampening)
numNewton = 10; % Number of Newton steps per iteration in parametric ss calculation


%% Initial Dynare run

% Set parameters and grids
setParameters;
computeGrids;
computePolynomials;

% Save all parameters
saveParameters;

delete steady_vars.mat;
rng(rng_seed);
dynare firstOrderDynamics_polynomials noclearall nopathchange;


%% Simulate data

if is_data_gen == 0
    
    % Load previous data
    load('simul.mat')
    load('simul_micro.mat');
    
else
    
    % Simulate macro variables
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data
    
    % Simulate micro variables
    sim_data_micro = simulate_micro(sim_struct, ts_micro, N_micro);
    save('simul_micro.mat', 'sim_data_micro');
    
end


%% Compute likelihood

loglikes = nan(length(param1_vals),length(param2_vals));
loglikes_macro = nan(length(param1_vals),length(param2_vals));
loglikes_micro = nan(length(param1_vals),length(param2_vals));

disp('Computing likelihood...');
timer_likelihood = tic;

if contains(pwd,'Laura')
    poolobj = parpool(2);
else
    poolobj = parpool;
end

for iter_i=1:length(param1_vals) % For each parameter...
    
    for iter_j=1:length(param2_vals) % For each parameter...
        
        % Set new parameters
        bbeta = param1_vals(iter_i);
        ttheta = param2_vals(iter_j);

        fprintf(['%s' repmat('%6.4f ',1,2),'%s\n'], '[bbeta,ttheta] = [',bbeta,ttheta,']');

        % Set parameters and grids
        setParameters;
        computeGrids;
        computePolynomials;

        saveParameters;         % Save all parameters
        setDynareParameters;    % Update Dynare parameters in model struct
        computeSteadyState;     % Compute steady state

        % Log likelihood of proposal
        [loglikes(iter_i,iter_j),loglikes_macro(iter_i,iter_j),loglikes_micro(iter_i,iter_j)] = ...
            loglike_compute('simul.mat', sim_data_micro, ts_micro, ...
                            num_smooth_draws, num_burnin_periods, ...
                            M_, oo_, options_);
    
    end
    
end

delete(poolobj);

likelihood_elapsed = toc(timer_likelihood);
fprintf('%s%8.2f\n', 'Done. Elapsed minutes: ', likelihood_elapsed/60);

cd('../../');
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');

% profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);