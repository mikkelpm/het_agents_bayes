function [sim_struct, sim_array] = simulate_model(num_sim_periods,num_burnin_periods,M_,oo_,options_)
    
    % Simulate from model
    options_.periods = num_sim_periods+num_burnin_periods;      % Total number of simulation periods
    endo_simul = simult(oo_.dr.ys, oo_.dr, M_, options_, oo_);  % Simulated data
    sim_array = endo_simul(:,num_burnin_periods+1:end)';        % Drop burn-in
    
    % Return structure with variable names
    sim_struct = struct;
    for j=1:M_.endo_nbr
        sim_struct.(deblank(M_.endo_names(j,:))) = sim_array(:,j);
    end

end