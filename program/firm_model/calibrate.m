% Calibrate heterogeneous firm model
% including settings for numerical model solution

% See Winberry (QE 2018)


%% Set economic parameters 

global ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
    bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ  cchi

% Technology
ttheta       = .256;                % capital coefficient
nnu          = .64;                 % labor coefficient
ddelta       = .085;                % depreciation (annual)
rrhoProd     = .53;                 % persistence of idiosyncratic shocks (annual)
ssigmaProd   = .0364;               % SD innovations of idiosycnratic shocks (annual)
aaUpper      = .011;                % no fixed cost region upper bound
aaLower      = -.011;               % no fixed cost region lower bound
ppsiCapital  = .0083;               % upper bound on fixed adjustment cost draws

% Preferences
bbeta        = .961;                % discount factor (annual)
ssigma       = 1;                   % coefficient of relative risk aversion
pphi         = 1 / 1e5;             % inverse Frisch elasticity of labor supply
nSS          = 1 / 3;               % hours worked in steady state
cchi         = 1;                   % labor disutility (will be calibrated in Dynare's steady state file to ensure labor supply = nSS)

% Aggregate shocks
rrhoTFP      = .859;                % persistence of aggregate TFP (annual)
ssigmaTFP    = .014;                % SD of innovations of aggregate TFP (annual)
rrhoQ        = .859;                % persistence of aggregate investment-specific shock (annual)
ssigmaQ      = .014;                % SD of innovations of aggregate investment-specific shock (annual)
corrTFPQ     = 0;                   % loading on TFP shock in evolution of investment-specific shock

%% Set approximation parameters

global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
    maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
    nCapitalQuadrature kRepSS wRepSS

% Order of approximation of value function
nProd        = 3;                   % order of polynomials in productivity
nCapital     = 5;                   % order of polynomials in capital
nState       = nProd * nCapital;    % total number of coefficients

% Shocks 
nShocks      = 3;                   % order of Gauss-Hermite quadrature over idiosyncratic shocks

% Finer grid for analyzing policy functions and computing histogram
nProdFine    = 60;
nCapitalFine = 40;
nStateFine   = nProdFine * nCapitalFine;

% Iteration on value function
maxIterations= 100;
tolerance    = 1e-6;
acc          = 500;                 % number of iterations in "Howard improvement step"
dampening    = 0;                   % weight on old iteration in updating step

% Approximation of distribution
nMeasure             = 2;           % order of polynomial approximating distribution
                                    % (NEED to change likelihood calculation if nMeasure ~= 2)
nProdQuadrature      = 8;           % number of quadrature points in productivity dimension
nCapitalQuadrature   = 10;          % number of quadrature points in capital dimension
nStateQuadrature     = nProdQuadrature * nCapitalQuadrature;
nMeasureCoefficients = (nMeasure * (nMeasure + 1)) / 2 + nMeasure;
