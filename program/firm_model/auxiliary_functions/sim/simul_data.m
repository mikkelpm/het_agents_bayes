% Simulate macro and micro data

sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data
    for j=1:nMeasureCoefficients
        sim_struct.(sprintf('%s%d', 'smpl_m', j)) = nan(T,1); % Set sample moments to missing everywhere
    end
    save_mat('simul', '-struct', 'sim_struct');                         % Save simulated data
    
    % draw micro data (incorporating truncation)
    simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, trunc_logn);
    save_mat('simul_data_micro','simul_data_micro');
    
    % Compute cross-sectional moments from micro data
    sim_struct_moments = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro);
    save_mat('simul_moments', '-struct', 'sim_struct_moments');
    
%     % draw individual productivities and incomes
%     simul_data_micro_indv_param = simulate_micro_indv_param(simul_data_micro);
%     save('simul_data_micro_indv_param.mat','simul_data_micro_indv_param');
