% Cov matrix of sample moments based on s.s.
smpl_m3 = oo_.mean(strcmp(oo_.var_list,'smpl_m3'));
smpl_m4 = oo_.mean(strcmp(oo_.var_list,'smpl_m4'));
smpl_m5 = oo_.mean(strcmp(oo_.var_list,'smpl_m5'));
N_micro_sample = mean(sum(~isnan(simul_data_micro(:,:,1)),2));
M_.H(3:end,3:end) = cov_smpl([smpl_m3 smpl_m4; smpl_m4 smpl_m5])/N_micro_sample;
