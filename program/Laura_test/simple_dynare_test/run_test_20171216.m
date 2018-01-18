clear all;


%% Initial Dynare run

dynare example1_new noclearall;         % Run Dynare once to process model file


%% Simulate data

% Settings
options_.periods = 100;                 % Simulate 100 periods
options_.drop = 0;                      % Don't drop any simulation periods

% Simulate
set_dynare_seed(20171215);              % Seed RNG
stoch_simul(var_list_);                 % Simulate

% Store simulated data
sim_vars = cellstr(options_.varobs)';   % Names of observed simulated
save('simul.mat', sim_vars{:});         % Save simulated data


%% Smoothing and likelihood

% Settings
options_.datafile = 'simul.mat';                % Data file
options_.periods = 0;                           % Don't simulate
options_.noprint = 1;                           % Don't print output after solving model
options_.mode_compute = 0;                      % Don't compute mode when "estimating"
options_.smoother = 1;                          % Return smoothed variables
options_.selected_variables_only = 1;           % Smooth selected variables (suppresses 
options_.smoothed_state_uncertainty = true;     % Return smoothing variance
options_.debug = true;                          % Returns likelihood value when running smoother

% Determine variables to smooth
smooth_vars = char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones

% Loop over alpha values
alphas = linspace(0.2, 0.5, 4); % alpha values to loop over
loglikes = zeros(size(alphas));

tic
for i=1:2%length(alphas) % For each alpha...

    fprintf('%s%4.2f\n', 'alpha=', alphas(i));
    set_param_value('alpha', alphas(i));                    % Set alpha
    disp('stoch_simul')
    stoch_simul(var_list_);                                 % Solve model with new parameters
    disp('estimation')
%     disp(['ss ' num2str(oo_.steady_state')])
    dynare_estimation(smooth_vars);                         % Run smoother
    loglikes(i) = -oo_.likelihood_at_initial_parameters;    % Store log likelihood value

end
% disp(['ss ' num2str(oo_.steady_state')])
t=toc