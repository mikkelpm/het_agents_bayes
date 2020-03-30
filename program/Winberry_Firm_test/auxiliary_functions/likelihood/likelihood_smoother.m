function ll = likelihood_smoother(dat_file, M_, oo_, options_)
  
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
    % NEED to better handle smooth_var
    [dataset_, dataset_info, xparam1, ~, M_new, options_new, oo_new, estim_params_,bayestopt_, bounds] = dynare_estimation_init('aggregateTFP', options_.dirname, 0, M_, options_, oo_, [], []);
    
    
    %% Likelihood
    
    % Compute likelihood (also re-solves model)
    disp('Computing likelihood...');
    [nll,~,~,~,~,~,~,~,~,~,~] = dsge_likelihood(xparam1,dataset_,dataset_info,options_new,M_new,estim_params_,bayestopt_,bounds,oo_new,[]);
    ll = -nll; % Log likelihood
   
end
