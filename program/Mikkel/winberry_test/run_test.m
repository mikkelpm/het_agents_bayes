clear all;


%% Settings

bbetas = linspace(0.9,0.99,10);     % beta values to loop over
T = 200;                            % Number of periods of simulated data
num_burnin_periods = 100;           % Number of burn-in periods for simulations
num_smooth_draws = 25;              % Number of draws from the smoothing distribution (for each beta)
rng_seed = 20180116;                % Random number generator seed for initial simulation


%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP;
	
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

% Aggregate Shocks
rrhoTFP = .859;										
ssigmaTFP = .014;


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

cd('./Auxiliary Functions');
saveParameters;


%% Initial Dynare run

dynare firstOrderDynamics_polynomials noclearall;                   % Run Dynare once to process model file


%% Simulate data

% Simulate
set_dynare_seed(rng_seed);                                          % Seed RNG
sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data

% Store simulated data
save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data


%% Smoothing and likelihood

% Determine variables to smooth
smooth_vars = ['w'; 'r']; %char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones

% Loop over beta values
loglikes = zeros(size(bbetas));
smooth_draws = cell(length(bbetas),num_smooth_draws);

for i=1:length(bbetas) % For each alpha...

    fprintf('%s%6.4f\n', 'beta=', bbetas(i));
    bbeta = bbetas(i);                          % Set beta
    saveParameters;
    setDynareParameters;

    % Durbin-Koopman simulation smoother
    [the_ll, ~, the_smooth_draws] = simulation_smoother('simul.mat', smooth_vars, num_smooth_draws, num_burnin_periods, M_, oo_, options_);
    loglikes(i) = the_ll;                       % Store log likelihood
    smooth_draws(i,:) = the_smooth_draws;       % Store smoothing draws

end

cd('../');
