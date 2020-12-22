% Simulate macro and micro data

sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate macro data
for i_Epsilon = 1:nEpsilon
    for i_Measure = 1:nMeasure
        sim_struct.(sprintf('%s%d%d', 'smpl_m', i_Epsilon, i_Measure)) = nan(T,1); 
            % Set sample moments to missing everywhere
    end
end
% Add measurement error to aggregate output
sim_struct.logAggregateOutput_noerror = sim_struct.logAggregateOutput;
sim_struct.logAggregateOutput = sim_struct.logAggregateOutput + ssigmaMeas*randn(T,1);
save_mat('simul', '-struct', 'sim_struct');                         % Save simulated data

% draw normalized individual incomes
simul_data_micro_aux = simulate_micro_aux(sim_struct, ts_micro, N_micro);

% draw individual productivities and incomes
simul_data_micro = simulate_micro(simul_data_micro_aux);
save_mat('simul_data_micro','simul_data_micro');

% Compute cross-sectional moments from micro data
sim_struct_moments = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro);
save_mat('simul_moments', '-struct', 'sim_struct_moments');

% Data set with missing third moment
sim_struct_moments2 = sim_struct_moments;
sim_struct_moments2.smpl_m13 = nan(T,1);
sim_struct_moments2.smpl_m23 = nan(T,1);
save_mat('simul_moments2', '-struct', 'sim_struct_moments2');

% Data set with missing second and third moments
sim_struct_moments1 = sim_struct_moments2;
sim_struct_moments1.smpl_m12 = nan(T,1);
sim_struct_moments1.smpl_m22 = nan(T,1);
save_mat('simul_moments1', '-struct', 'sim_struct_moments1');
