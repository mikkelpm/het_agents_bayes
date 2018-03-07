function simul_data_hh_indv_param = simulate_hh_indv_param(simul_data_hh)

% Simulate household employment and income
global mu_l;
disp('Simulating household employment and income (individual productivity)...');

simul_data_hh_indv_param = simul_data_hh;
simul_data_hh_indv_param(:,:,2) ...
    = lognrnd(mu_l,sqrt(-2*mu_l),size(simul_data_hh_indv_param(:,:,2)))...
    .*simul_data_hh_indv_param(:,:,2); % here sample from lognomral directly

end