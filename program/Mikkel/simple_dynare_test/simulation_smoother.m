function [ll, smooth_means, smooth_draws, oo_new] = simulation_smoother(dat_file, smooth_vars, num_smooth_draws, num_burnin_periods, M_, oo_, options_)
  
    % Likelihood calculation and Durbin-Koopman simulation smoother


    %% Preliminaries
    
    % Set options
    options_.datafile = dat_file;                   % Data file
    options_.smoother = 1;                          % Return smoothed variables
    options_.noprint = 1;                           % Don't print output when solving model
    options_.mode_compute = 0;                      % Don't compute mode when "estimating"
    options_.selected_variables_only = 1;           % Smooth selected variables (suppresses a prompt)
    
    % Parameter, model, and data information
    [dataset_, dataset_info, xparam1, ~, M_, options_, oo_, estim_params_,bayestopt_, bounds] = dynare_estimation_init(smooth_vars, options_.dirname, 0, M_, options_, oo_, [], []);
    
    
    %% Likelihood
    
    % Compute likelihood (also re-solves model)
    [nll,~,~,~,~,~,~,~,~,~,oo_new] = dsge_likelihood(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,bounds,oo_,[]);
    ll = -nll; % Log likelihood
    
    
    %% Conditional mean smoothing
    
    % Run smoother (also re-solves model)
    oo_smooth = mean_smoother(dataset_.data',xparam1,dataset_,dataset_info,M_,oo_new,options_,bayestopt_,estim_params_);
    smooth_means = oo_smooth.SmoothedVariables; % Smoothed means
    
    
    %% Simulation smoother
    
    T = dataset_.nobs;  % Sample size
    smooth_draws = cell(1,num_smooth_draws);
    
    for iter_i=1:num_smooth_draws % For each draw...
        
        % Simulate fake data
        [sim_results,sim_array] = simulate_model(T,num_burnin_periods,M_,oo_new,options_);
        
        % Run mean smoother on fake data
        sim_oo_smooth = mean_smoother(sim_array(:,options_.varobs_id)',xparam1,dataset_,dataset_info,M_,oo_new,options_,bayestopt_,estim_params_);
        
        % Create draws from smoothing distribution
        % See Durbin & Koopman, 2012, 2nd ed, ch. 4.9.1
        smooth_draws{iter_i} = struct;
        for iter_j=1:size(smooth_vars,1)
            the_var = deblank(smooth_vars(iter_j,:));
            smooth_draws{iter_i}.(the_var) = smooth_means.(the_var) + sim_results.(the_var) - sim_oo_smooth.SmoothedVariables.(the_var);
        end

    end

end