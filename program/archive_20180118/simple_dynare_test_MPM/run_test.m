clear all;


%% Settings

alphas = linspace(0.2, 0.5, 10);    % alpha values to loop over
T = 200;                            % Number of periods of simulated data
num_burnin_periods = 100;           % Number of burn-in periods for simulations
num_smooth_draws = 20;              % Number of draws from the smoothing distribution (for each alpha)
rng_seed = 20171219;                % Random number generator seed


%% Initial Dynare run

dynare example1_new noclearall;         % Run Dynare once to process model file


%% Simulate data

% Simulate
set_dynare_seed(rng_seed);                                          % Seed RNG
sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data

% Store simulated data
save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data


%% Smoothing and likelihood

% Determine variables to smooth
smooth_vars = char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones

% Loop over alpha values
loglikes = zeros(size(alphas));
smooth_draws = cell(length(alphas),num_smooth_draws);

for i=1:length(alphas) % For each alpha...

    fprintf('%s%4.2f\n', 'alpha=', alphas(i));
    set_param_value('alpha', alphas(i));        % Set alpha

    % Durbin-Koopman simulation smoother
    [the_ll, ~, the_smooth_draws, the_oo_new] = simulation_smoother('simul.mat', smooth_vars, num_smooth_draws, num_burnin_periods, M_, oo_, options_);
    loglikes(i) = the_ll;                       % Store log likelihood
    smooth_draws(i,:) = the_smooth_draws;       % Store smoothing draws
    
    % Update results with new steady state
    oo_ = the_oo_new;

end