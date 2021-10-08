% Simulate macro and micro data

sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate macro data
save_mat('simul', '-struct', 'sim_struct');                         % Save simulated data

% draw micro data (incorporating truncation)
simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, trunc_logn);
save_mat('simul_data_micro','simul_data_micro');
