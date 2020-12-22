function [prop_draw,is_adapt] = rwmh_propose(curr_draw, stepsize, chol, p_adapt, s_diffused)

% Random Walk Metropolis-Hastings proposal draw

is_adapt = rand(1)<p_adapt;
if is_adapt
    prop_draw = curr_draw + stepsize*randn(size(curr_draw))*chol; % Adaptive component
else
    prop_draw = curr_draw + s_diffused*randn(size(curr_draw)); % Diffused component
end

end