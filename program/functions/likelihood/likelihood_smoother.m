function [ll, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] = likelihood_smoother(dat_file, smooth_vars, M_, oo_, options_, do_smooth)
  
    % Likelihood and mean smoother


    %% Preliminaries
    
    % Set options
    options_.datafile = dat_file;                   % Data file
    options_.smoother = 1;                          % Return smoothed variables
    options_.noprint = 1;                           % Don't print output when solving model
    options_.debug = 0;
    options_.mode_compute = 0;                      % Don't compute mode when "estimating"
    options_.selected_variables_only = 1;           % Smooth selected variables (suppresses a prompt)
    
    % Parameter, model, and data information 
    [dataset_, dataset_info, xparam1, ~, M_new, options_new, oo_new, estim_params_,bayestopt_, bounds] = dynare_estimation_init(smooth_vars, options_.dirname, 0, M_, options_, oo_, [], []);
    
    
    %% Likelihood
    
    % Compute likelihood (also re-solves model)
    disp('Computing likelihood...');
        [nll,~,~,~,~,~,~,~,~,~,oo_new] = dsge_likelihood(xparam1,dataset_,dataset_info,options_new,M_new,estim_params_,bayestopt_,bounds,oo_new,[]);
    ll = -nll; % Log likelihood
    
    
    %% Conditional mean smoothing
    
    smooth_means = [];
    if do_smooth
        % Run smoother
        disp('Mean smoother...');
        oo_smooth = mean_smoother(dataset_.data',xparam1,dataset_,dataset_info,M_new,oo_new,options_new,bayestopt_,estim_params_);
        smooth_means = oo_smooth.SmoothedVariables;
    end
   
end
