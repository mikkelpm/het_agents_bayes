function [sim_struct, sim_array] = simulate_model(num_sim_periods,num_burnin_periods,M_,oo_,options_,varargin)
    
    % Simulate macro data from Dynare model

    % Check if shocks are already supplied
    if ~isempty(varargin) && ~isempty(varargin{1})
        shocks = varargin{1};
    else
        shocks = simulate_shocks(M_, num_sim_periods+num_burnin_periods, 1);
    end

    % Simulate from model
    endo_simul = simult_(M_, options_, oo_.dr.ys, oo_.dr, shocks, 1);   % Simulated data
    sim_array = endo_simul(:,num_burnin_periods+2:end)';                % Drop burn-in (and initial condition)
    
    % Return structure with variable names
    sim_struct = struct;
    for j=1:M_.endo_nbr
        sim_struct.(deblank(M_.endo_names{j})) = sim_array(:,j);
    end

end