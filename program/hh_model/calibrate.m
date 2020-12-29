% Calibrate heterogeneous household model
% including settings for numerical model solution

% See Winberry (QE 2018)


%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
    mmu rrhoTFP ssigmaTFP ttau mu_l ssigmaMeas;
  
% Preferences
bbeta = .96;                    % discount factor (annual calibration)
ssigma = 1;                     % coefficient of relative risk aversion
aaBar = 0;                      % borrowing constraint

% Technology
aalpha = .36;                   % capital share
ddelta = .1;                    % depreciation rate (annual calibration)

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
splineOpt = 0; % MUST BE SET TO 0 (for polynomial approximation)

% Whether to print out results from steady state computation
displayOpt = 'off'; % 'iter-detailed' or 'off'

% Order of approximation
nEpsilon = 2;
nAssets = 25; % number of gridpoints in spline approximation or polynomials in polynomial approximation

% Finer grid for analyzing policy functions
nAssetsFine = 100;

% Approximation of distribution
nMeasure = 3; % NEED to change likelihood calculation if nMeasure ~= 3
nAssetsQuadrature = 8;

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;
