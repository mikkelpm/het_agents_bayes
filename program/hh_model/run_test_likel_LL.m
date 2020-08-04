clear all;
addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');

%% Settings

% Decide what to do
is_data_gen = 1; % whether simulate data:  
                 % 0: no simulation
                 % 1: simulation

% Model/data settings
T = 100;                                % Number of periods of simulated macro data
ts_micro = 10:10:T;                        % Time periods where we observe micro data
N_micro = 1e2;                             % Number of households per non-missing time period
micro_mom_order = 2;                    % Observe micro moments up to this order for inference

% Parameter values to check
l_param = {'bbeta','ssigmaMeas','mu_l'};
tex_param = {'\beta','\sigma_e','\mu_\lambda'};
n_param = length(l_param);
param_vals_mult = unique([1 linspace(0.75,1.25,51)]); % Multiples of true parameter to compute

% Likelihood settings
num_smooth_draws = 500;                 % Number of draws from the smoothing distribution (for unbiased likelihood estimate)
num_interp = 100;                       % Number of interpolation grid points for calculating density integral

% Numerical settings
num_burnin_periods = 100;               % Number of burn-in periods for simulations
rng_seed = 20180727;                    % Random number generator seed for initial simulation
rng('default');


%% Set economic parameters

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
    mmu rrhoTFP ssigmaTFP ttau mu_l ssigmaMeas;

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
splineOpt = 0;	% MUST BE SET TO 0

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

cd('./auxiliary_functions/dynare');

delete steady_vars.mat;
saveParameters;
    
    
%% Initial Dynare run

rng(rng_seed);
dynare firstOrderDynamics_polynomials noclearall nopathchange; % Run Dynare once to process model file


%% Simulate data

if is_data_gen == 1

    % Simulate
    set_dynare_seed(rng_seed);                                          % Seed RNG
    sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    for i_Epsilon = 1:nEpsilon
        for i_Measure = 1:nMeasure
            sim_struct.(sprintf('%s%d%d', 'smpl_m', i_Epsilon, i_Measure)) = nan(T,1); 
                % Set sample moments to missing everywhere
        end
    end
    % Add measurement error to aggregate output
    sim_struct.logAggregateOutput_noerror = sim_struct.logAggregateOutput;
    sim_struct.logAggregateOutput = sim_struct.logAggregateOutput + ssigmaMeas*randn(T,1);
    save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data

    % draw normalized individual incomes
    simul_data_micro_aux = simulate_micro_aux(sim_struct, ts_micro, N_micro);

    % draw individual productivities and incomes
    simul_data_micro = simulate_micro(simul_data_micro_aux);
    save('simul_data_micro.mat','simul_data_micro');

    % Compute cross-sectional moments from micro data
    sim_struct_moments = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro);
    save('simul_moments.mat', '-struct', 'sim_struct_moments');
    
    for im=micro_mom_order+1:nMeasure
        % Set higher-order sample moments to missing
        sim_struct_moments.(sprintf('%s%d', 'smpl_m1', im)) = nan(T,1);
        sim_struct_moments.(sprintf('%s%d', 'smpl_m2', im)) = nan(T,1);
    end
    save('simul_moments2.mat', '-struct', 'sim_struct_moments');

else

    % Load previous data
%         load('simul.mat')
    load('simul_data_micro.mat');

end


%% Loop over parameters to estimate

if isempty(gcp)
    parpool;
end

for i_param = 1:n_param
    
    param_plot = l_param{i_param}; % Parameter to vary when plotting likelihood
    
    
    %% Compute likelihood
    
    eval(['param_true = ' param_plot ';']);
    param_vals = param_true*param_vals_mult;
    loglikes = nan(length(param_vals),3);
    loglikes_macro = nan(length(param_vals),3);
    loglikes_micro = nan(length(param_vals),3);
    
    disp('Computing likelihood...');
    timer_likelihood = tic;
    
    for iter_i=1:length(param_vals) % For each parameter value...
            
        % Set new parameter
        eval([param_plot ' = param_vals(iter_i);']);

        fprintf('%d%s%d%s%s%s%6.4f\n', iter_i, '/', length(param_vals), ': ', param_plot, ' = ', eval(param_plot));

        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        try
            
            compute_steady_state;   % Compute steady state
            compute_meas_err;       % Update measurement error var-cov matrix
        
        catch ME
            
            disp('Error encountered. Message:');
            disp(ME.message);
            continue
            
        end
        
        for i_setup = 1:3
            try
                if i_setup == 1
                    % Macro + full info micro
                    [loglikes(iter_i,1), loglikes_macro(iter_i,1), loglikes_micro(iter_i,1)] = ...
                        loglike_compute('simul.mat', simul_data_micro, ts_micro, ...
                        num_smooth_draws, num_interp, num_burnin_periods, ...
                        M_, oo_, options_);
                elseif i_setup == 2
                    % Macro + moments w/ SS meas. err.
                    [loglikes(iter_i,2), loglikes_macro(iter_i,2), loglikes_micro(iter_i,2)] = ...
                        loglike_compute('simul_moments.mat', [], ts_micro, ...
                        0, num_interp, num_burnin_periods, ...
                        M_, oo_, options_);
                elseif i_setup == 3
                    % Macro + moments w/ SS meas. err. (observe up to 2nd
                    % moment)
                    [loglikes(iter_i,3), loglikes_macro(iter_i,3), loglikes_micro(iter_i,3)] = ...
                        loglike_compute('simul_moments2.mat', [], ts_micro, ...
                        0, num_interp, num_burnin_periods, ...
                        M_, oo_, options_);
%                 else
%                     % Macro + moments w/o meas. err.
%                     M_.H(2:end,2:end) = 1e-8*eye(nMeasure*nEpsilon); % Only a little bit of meas. err. to avoid singularity
%                     [loglikes(iter_i,3), loglikes_macro(iter_i,3), loglikes_micro(iter_i,3)] = ...
%                         loglike_compute('simul_moments.mat', [], ts_micro, ...
%                         0, num_interp, num_burnin_periods, ...
%                         M_, oo_, options_);
                end
            catch ME
                
                disp('Error encountered. Message:');
                disp(ME.message);
                
            end
        end
        
    end
    
    likelihood_elapsed = toc(timer_likelihood);
    fprintf('%s%8.2f\n', 'Done. Elapsed minutes: ', likelihood_elapsed/60);
    
    %% Plot

    loglikes_plot = real([loglikes loglikes_macro(:,1)]);
    plot(param_vals,loglikes_plot-max(loglikes_plot,[],1),'o-');
    the_ylim = ylim;
    line(param_true*ones(1,2), the_ylim, 'Color', 'k', 'LineStyle', ':');
    ylim(the_ylim);
    legend({'FI micro', '3rd moments micro', '2nd moment micro', 'macro'}, 'Location', 'SouthWest');
    title(tex_param{i_param})
    savefig(l_param{i_param})
    
    %% Change back to true parameter value
    
    eval([param_plot ' = param_true;']);
    
    %% Save
    
    save([l_param{i_param} '.mat'],'loglikes','loglikes_micro','loglikes_macro','param_vals')
    
end

delete(gcp('nocreate'));

cd('../../');
rmpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');