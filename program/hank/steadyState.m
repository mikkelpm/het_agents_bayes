% Computes and analyzes steady state with no aggregate shocks
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

% profile on
oldFolder = cd('./auxiliary_functions');


%% Define parameters

% ECONOMIC PARAMETERS

% Preferences
bbeta = .96; % Discount factor (annual calibration)
ppsi = 1;

% Taxes
ttau = 0.1;
mmu = 0.05;

% Idioynscratic productivity
vzGrid = [0;1];
mzTransition = [0.5 0.5; 0.5 0.5];

% Aggregate productivity
AAss = 1; % Steady state productivity level

% APPROXIMATION PARAMETERS

% Whether to print out results from steady state computation
displayOpt = 'off';       % 'iter-detailed' or 'off'

% Order of approximation
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


%% Set remaining parameters

setParameters;


%% Compute Steady State

% Solve for steady state
coreSteadyState;
profile viewer;
% profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);


cd(oldFolder);