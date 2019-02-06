% Computes and analyzes steady state with no aggregate shocks
% 2019-02-04

clear all;
close all;

% profile on
oldFolder = cd('./auxiliary_functions');


%% Define parameters

% ECONOMIC PARAMETERS

% Preferences
bbeta = .96^(1/4); % Discount factor (quarterly)
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

% Idioynscratic productivity
vzGrid = [0.5;1.5];
mzTransition = [0.7 0.3; 0.1 0.9];

% Aggregate productivity
A_SS = 3; % Steady state productivity level

% APPROXIMATION PARAMETERS

% Whether to print out results from steady state computation
displayOpt = 'iter-detailed';       % 'iter-detailed' or 'off'

% Order of approximation
nAssets = 25; % number of polynomials in polynomial approximation

% Finer grid for analyzing policy functions
nAssetsFine = 100;

% Approximation of distribution
nMeasure = 3;
nAssetsQuadrature = 8;

% Steady state
tolerance_SS_root = 0.01; % Numerical tolerance for root finding
tolerance_SS_invhist = 1e-12; % Numerical tolerance for invariant distribution of histogram approach

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;


%% Set remaining parameters

setParameters;


%% Compute Steady State

% Solve for steady state
coreSteadyState;
% profile viewer;
% profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);


cd(oldFolder);