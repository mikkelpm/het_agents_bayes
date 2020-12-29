% Random Walk Metropolis-Hastings (RWMH) MCMC procedure

accepts = zeros(mcmc_num_iter,1);

mcmc_num_draws = floor(mcmc_num_iter/mcmc_thin); % Number of stored draws
post_draws = nan(mcmc_num_draws,length(param_names));

curr_draw = mcmc_init;
curr_logpost = -Inf;
loglikes_prop = nan(mcmc_num_draws,1);
loglikes_prop_macro = nan(mcmc_num_draws,1);
loglikes_prop_micro = nan(mcmc_num_draws,1);

the_stepsize = mcmc_stepsize_init;      % Initial RWMH step size
the_stepsize_iter = 1;
the_chol = eye(length(curr_draw));      % Initial RWMH proposal var-cov matrix

fprintf('\nMCMC...\n');
timer_mcmc = tic;

for i_mcmc=1:mcmc_num_iter % For each MCMC step...
    
    print_param(transf_to_param(curr_draw), param_names, 'current');
    
    % Proposed draw (modified to always start with initial draw)
    [prop_draw, is_adapt] = rwmh_propose(curr_draw, (i_mcmc>1)*the_stepsize, the_chol, mcmc_p_adapt, mcmc_stepsize_init); % Proposal
    update_param(transf_to_param(prop_draw), param_names);
    print_param(transf_to_param(prop_draw), param_names, 'proposed');
    
    try
        
        [the_loglike_prop, the_loglike_prop_macro, the_loglike_prop_micro] = ...
            ll_fct(M_, oo_, options_);
        
        % Log prior density of proposal
        logprior_prop = prior_logdens_transf(prop_draw);

        % Accept/reject
        [curr_draw, curr_logpost, accepts(i_mcmc), the_log_ar] = rwmh_accrej(curr_draw, prop_draw, curr_logpost, logprior_prop+the_loglike_prop);
    
    catch ME
        
        disp('Error encountered. Message:');
        disp(ME.message);
        
        the_loglike_prop = nan;
        the_loglike_prop_macro = nan;
        the_loglike_prop_micro = nan;
        
    end
    
    % Store (with thinning)
    if mod(i_mcmc, mcmc_thin)==0
        post_draws(i_mcmc/mcmc_thin,:) = transf_to_param(curr_draw);
        loglikes_prop(i_mcmc/mcmc_thin) = the_loglike_prop;
        loglikes_prop_macro(i_mcmc/mcmc_thin) = the_loglike_prop_macro;
        loglikes_prop_micro(i_mcmc/mcmc_thin) = the_loglike_prop_micro;
    end
    
    % Print acceptance rate
    fprintf('%s%5.1f%s\n', 'Accept. rate last 100: ', 100*mean(accepts(max(i_mcmc-99,1):i_mcmc)), '%');
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_mcmc, '/', mcmc_num_iter);
    
    % Adapt proposal step size
    if is_adapt
        if exist('the_log_ar','var')
            the_ar = exp(the_log_ar);
        else
            the_ar = 0;
        end
        [the_stepsize, the_stepsize_iter] = adapt_stepsize(the_stepsize, the_stepsize_iter, i_mcmc, the_ar, mcmc_c, mcmc_ar_tg);
    end
        
    % Adapt proposal covariance matrix
    [the_chol, the_stepsize_iter] = adapt_cov(the_chol, the_stepsize_iter, mcmc_adapt_iter, i_mcmc, post_draws, mcmc_thin, mcmc_adapt_diag, mcmc_adapt_param);
    
    % Save middle steps in case reach time limit on the server
    if mod(i_mcmc,1000) == 0 && i_mcmc < mcmc_num_iter
        save_mat(fullfile(save_folder, model_name));
    end
    
end

rng_status = rng;
mcmc_elapsed = toc(timer_mcmc);
fprintf('%s%8.2f\n', 'MCMC done. Elapsed minutes: ', mcmc_elapsed/60);
