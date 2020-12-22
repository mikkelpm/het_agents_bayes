function smooth_draw = simulation_smoother(smooth_means, smooth_vars, num_burnin_periods, sim_shocks, M_, oo_, options_, dataset_, dataset_info, xparam1, estim_params_, bayestopt_)
  
    % Simulation smoother
    % Durbin & Koopman (Biometrika 2002)
    
    
    T = dataset_.nobs;  % Sample size
    
    % Simulate fake data
    [sim_results,sim_array] = simulate_model(T,num_burnin_periods,M_,oo_,options_,sim_shocks);

    % Run mean smoother on fake data
    sim_oo_smooth = mean_smoother(sim_array(:,options_.varobs_id)',xparam1,dataset_,dataset_info,M_,oo_,options_,bayestopt_,estim_params_);
    sim_smooth_means = sim_oo_smooth.SmoothedVariables;

    % Create draws from smoothing distribution
    % See Durbin & Koopman, 2012 book, 2nd ed, ch. 4.9.1
    smooth_draw = struct;
    for iter_j=1:size(smooth_vars,1)
        the_var = deblank(smooth_vars{iter_j});
        smooth_draw.(the_var) = smooth_means.(the_var) + sim_results.(the_var) - sim_smooth_means.(the_var);
    end


end
