% Cov matrix of sample moments based on s.s.
smpl_m12 = oo_.mean(strcmp(oo_.var_list,'smpl_m12'));
smpl_m22 = oo_.mean(strcmp(oo_.var_list,'smpl_m22'));
M_.H(1,1) = ssigmaMeas^2;
M_.H(2:end,2:end) = cov_smpl([smpl_m12 smpl_m22])/N_micro;