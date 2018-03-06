function simul_data_hh_indv_param = simulate_hh_indv_param(simul_data_hh)

% Simulate household employment and income
global fl_param;
disp('Simulating household employment and income (individual productivity)...');

simul_data_hh_indv_param = simul_data_hh;
simul_data_hh_indv_param(:,:,2) ...
    = exp(lognrnd(fl_param.mm(1),sqrt(-fl_param.gg(2)/0.5),size(simul_data_hh_indv_param(:,:,2))))...
    .*simul_data_hh_indv_param(:,:,2); % here sample from lognomral directly

end