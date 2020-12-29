function sim_struct = simulate_micro_moments(sim_struct, simul_data_micro, T, ts_micro)

    % Compute cross-sectional moments
    
    the_simul = nan(T,size(simul_data_micro,2),2);
    the_simul(ts_micro,:,1) = simul_data_micro(:,:,1); % Fill in non-missing obs.
    the_simul(ts_micro,:,2) = simul_data_micro(:,:,2);
    
    % Add moments to simulation struct
    ix0 = the_simul(:,:,1)==0;
    vN0 = sum(ix0,2);
    ix1 = the_simul(:,:,1)==1;
    vN1 = sum(ix1,2);

    % For moment-based methods up to the 3rd moments. 
    % NEED to change for moment-based methods with higher order moments.
    sim_struct.smpl_m11 = sum(the_simul(:,:,2).*ix0,2)./vN0;
    sim_struct.smpl_m12 = sum((the_simul(:,:,2)-sim_struct.smpl_m11).^2.*ix0,2)./vN0;
    sim_struct.smpl_m13 = sum((the_simul(:,:,2)-sim_struct.smpl_m11).^3.*ix0,2)./vN0;
    sim_struct.smpl_m21 = sum(the_simul(:,:,2).*ix1,2)./vN1;
    sim_struct.smpl_m22 = sum((the_simul(:,:,2)-sim_struct.smpl_m21).^2.*ix1,2)./vN1;
    sim_struct.smpl_m23 = sum((the_simul(:,:,2)-sim_struct.smpl_m21).^3.*ix1,2)./vN1;

end