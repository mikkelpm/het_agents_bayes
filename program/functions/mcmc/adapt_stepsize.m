function [new_stepsize, new_iter] = adapt_stepsize(curr_stepsize, curr_stepsize_iter, iter, curr_ar, c, ar_tg)

    % Adaptive step size, cf. Atchade and Rosenthal (2005), Section 4.1
    new_stepsize = exp(log(curr_stepsize)+(iter>1)*curr_stepsize_iter^(-c)*(curr_ar-ar_tg)); 
    new_stepsize = min(max(new_stepsize,1e-10),1e10);
    new_iter = curr_stepsize_iter + 1;
    fprintf('%s%8.6f\n', 'New step size: ', new_stepsize);
        
end