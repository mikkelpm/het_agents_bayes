function sim_struct = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro)

    % Compute cross-sectional moments
    
    the_simul = nan(T,size(simul_data_micro,2),size(simul_data_micro,3));
    the_simul(ts_micro,:,:) = simul_data_micro; % Fill in non-missing obs.
    
    % Add moments to simulation struct
    sim_struct.smpl_m1 = nanmean(the_simul(:,:,1),2);
    sim_struct.smpl_m2 = nanmean(the_simul(:,:,2),2);
    sim_struct.smpl_m3 = nanvar(the_simul(:,:,1),1,2);
    sim_struct.smpl_m4 = nanmean((the_simul(:,:,1)-sim_struct.smpl_m1).*(the_simul(:,:,2)-sim_struct.smpl_m2),2);
    sim_struct.smpl_m5 = nanvar(the_simul(:,:,2),1,2);

end