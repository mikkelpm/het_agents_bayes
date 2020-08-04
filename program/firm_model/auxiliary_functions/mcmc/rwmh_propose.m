function prop_draw = rwmh_propose(curr_draw, stepsize, chol)

    % Proposed draw
    prop_draw = curr_draw + stepsize*randn(size(curr_draw))*chol; % Proposal

end