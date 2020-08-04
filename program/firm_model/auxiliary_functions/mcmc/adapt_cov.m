function [new_chol, new_stepsize_iter] = adapt_cov(curr_chol, curr_stepsize_iter, adapt_iter, iter, post_draws, thin, adapt_diag, adapt_param)

    % Update RWMH proposal var-cov matrix

    new_chol = curr_chol;
    new_stepsize_iter = curr_stepsize_iter;

    the_indx = find(adapt_iter==iter,1); % Find current iteration index in list of adaptation iterations
    if ~isempty(the_indx) && the_indx>1 % If in list, but not first...
        the_start = 1 + (the_indx>1)*adapt_iter(the_indx-1)/thin; % Start of current adaptation window
        the_end = iter/thin; % End of current adaptation window
        the_cov = cov(post_draws(the_start:the_end,:)); % Var-cov matrix over current adaptation window
        the_cov = the_cov/(mean(sqrt(diag(the_cov)))^2); % Normalize average std deviation
        the_n = the_end-the_start+1;
        the_cov_shrink = (the_n*the_cov + adapt_param*eye(size(post_draws,2))) ...
                         /(the_n+adapt_param); % Shrink var-cov matrix toward identity
        if adapt_diag % If adapt only to std devs...
            new_chol = diag(sqrt(diag(the_cov_shrink)));
        else % If adapt full var-cov matrix...
            new_chol = chol(the_cov_shrink);
        end
        disp('New RWMH proposal var-cov matrix:');
        disp(new_chol'*new_chol);
        disp('Square root of diagonal');
        disp(sqrt(diag(new_chol'*new_chol)));
        new_stepsize_iter = 1; % Reset stepsize adaptation
    end

end