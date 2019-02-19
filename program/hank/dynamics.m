% Shell which declares parameters and calls Dynare to solve model using
% first order approximation of aggregate dynamics
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

% profile on

cd('./auxiliary_functions');

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
% vzGrid = [0.5;1.5];
% mzTransition = [0.8 0.2;
%                 0.2 0.8]; % (i,j) element: P(z'=z_j | z=z_i)

% Aggregate productivity
A_SS = 3; % Steady state aggr productivity level

% Equity shares
% vShareGrid = [0; 1; 2]; % Profit shares for each household type
% vShareFraction = [1/3; 1/3; 1/3]; % Fractions of each household type
vShareGrid = 1; %[0.5; 1.5]; % Profit shares for each household type
vShareFraction = 1; %[1/2; 1/2]; % Fractions of each household type

% Aggregate Shocks
rrhoTFP = .95; % Autocorrelation of TFP
ssigmaTFP = .005; % Std dev of innovation to TFP shock
rrhoMP = 0.9; % Autocorrelation of monetary policy shock
ssigmaMP = 0.0025; % Std dev of innovation to monetary policy shock

% Dynamic Parameters (KMV (2018))
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
nMeasure = 3;
nAssetsQuadrature = 8;

% Steady state
tolerance_SS_root = 0.001; % Numerical tolerance for root finding
tolerance_SS_invhist = 1e-12; % Numerical tolerance for invariant distribution of histogram approach
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .5;%.95;
numNewton = 10; % Number of Newton steps per iteration in parametric ss calculation

%%
%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (use polynomials to approximate conditional expectation)
computePolynomials;

%----------------------------------------------------------------
% Save parameters in .mat files to import into Dynare 
%----------------------------------------------------------------

% Economic parameters
save economicParameters.mat bbeta ppsi nnu bbBar eepsilon ttau vvarthetaB vvarthetaT...
    vzGrid vzInvariant mzTransition A_SS w_SS r_RepSS N_RepSS vShareGrid vShareFraction...
    rrhoTFP ssigmaTFP rrhoMP ssigmaMP ttheta pphi

% Approximation parameters
save approximationParameters.mat nz nAssets nState nAssetsFine nStateFine...
    nAssetsQuadrature nStateQuadrature nMeasure nMeasureCoefficients nShare...
    assetsMin assetsMax...
    tolerance_SS_root tolerance_SS_invhist maxIterations tolerance dampening numNewton

% Grids
save grids.mat vAssetsGridZeros vAssetsGrid mzGrid mAssetsGrid mzPrimeGrid vAssetsGridFine ...
    vAssetsGridFineZeros mzGridFine mAssetsGridFine mzPrimeGridFine vQuadratureWeights ...
    vAssetsGridQuadratureZeros vAssetsGridQuadrature mzGridQuadrature mAssetsGridQuadrature

% Polynomials
save polynomials.mat vAssetsPoly vAssetsPolySquared vAssetsPolyFine vAssetsPolyQuadrature vAssetsPolyBC

%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

dynare firstOrderDynamics_polynomials noclearall nopathchange;

cd('../')

% profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);