% Part of measurement error cov matrix that doesn't change over parameter values

% Here: cov matrix of cross-sectional moments (assumed constant over time)
% estimated using sample higher moments

smpl_m = zeros(2,6);
num_obs = zeros(2,1);
for eepsilon=0:1 % For each employment status
    aux = simul_data_micro(:,:,2);
    aux(simul_data_micro(:,:,1)~=eepsilon) = nan; % Keep only data for current employment status
    smpl_m(eepsilon+1,1) = nanmean(aux(:)); % Grand mean
    aux_demean = aux-nanmean(aux,1); % Subtract time-specific means
    smpl_m(eepsilon+1,2:6) = nanmean(aux_demean(:).^(2:6)); % Other moments
    num_obs(eepsilon+1) = mean(sum(simul_data_micro(:,:,1)==eepsilon,2)); % Average number of observations in each period
end
clearvars aux aux_demean;

M_.H(2:4,2:4) = cov_smpl(smpl_m(1,:))/num_obs(1);
M_.H(5:7,5:7) = cov_smpl(smpl_m(2,:))/num_obs(2);