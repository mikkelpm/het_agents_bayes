function simul_data_micro = simulate_micro(simul_data_micro_aux)

% Simulate household employment and income
global mu_l;
disp('Simulating household employment and income (individual productivity)...');

simul_data_micro = simul_data_micro_aux(:,:,[1 2]);
simul_data_micro(:,:,2) ...
    = lognrnd(mu_l,sqrt(-2*mu_l),size(simul_data_micro(:,:,2)))...
    .*simul_data_micro(:,:,2); % here sample from lognormal directly

end