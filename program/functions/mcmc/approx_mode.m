% Grid search for finding approximate posterior mode

optim_numgrid = size(optim_grid,1);
optim_logpost = nan(optim_numgrid,1);

fprintf('\nOptimization...\n');
timer_optim = tic;

for i_optim=1:optim_numgrid % Cycle through grid points

    the_param = optim_grid(i_optim,:);
    print_param(the_param, param_names, 'current');
    update_param(the_param, param_names);

    try
        the_loglike = ll_fct(M_, oo_, options_);
        the_logprior = prior_logdens_transf(param_to_transf(the_param));
        optim_logpost(i_optim) = the_loglike + the_logprior;
    catch ME
        disp('Error encountered. Message:');
        disp(ME.message);
    end

    % Print progress
    fprintf('%s%6d%s%6d\n\n', 'Progress: ', i_optim, '/', optim_numgrid);

end

% Find maximum of log posterior on grid
[~,optim_maxind] = max(optim_logpost);
optim_maxparam = optim_grid(optim_maxind,:);

optim_elapsed = toc(timer_optim);
fprintf('%s%8.2f\n', 'Optimization done. Elapsed minutes: ', optim_elapsed/60);

print_param(optim_maxparam, param_names, 'approximate mode');

% Start MCMC at approximate mode
mcmc_init = param_to_transf(optim_maxparam);
