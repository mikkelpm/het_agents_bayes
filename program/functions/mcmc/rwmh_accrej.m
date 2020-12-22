function [new_draw, new_logpost, accept, log_ar] = rwmh_accrej(curr_draw, prop_draw, curr_logpost, prop_logpost)

    % Accept/reject Metropolis-Hastings proposal

    log_ar = min(prop_logpost-curr_logpost,0);
    
    if log_ar > log(rand())
        fprintf('%s\n', 'Accepted!');
        new_draw = prop_draw;
        new_logpost = prop_logpost;
        accept = 1;
    else
        fprintf('%s\n', 'Rejected!');
        new_draw = curr_draw;
        new_logpost = curr_logpost;
        accept = 0;
    end

end