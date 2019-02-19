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
vShareGrid = [0; 1; 2]; % Profit shares for each household type
vShareFraction = [1/3; 1/3; 1/3]; % Fractions of each household type
% vShareGrid = 1;
% vShareFraction = 1;

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


%% Set parameters

% Set additional parameters
setParameters;

% Grids
computeGrids;

% Polynomials over grids (use polynomials to approximate conditional expectation)
computePolynomials;

% Save all parameters
saveParameters;


%% Run Dynare

dynare firstOrderDynamics_polynomials noclearall nopathchange;

cd('../../');
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');

% profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);