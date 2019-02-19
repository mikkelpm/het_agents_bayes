function [loglike, loglike_macro, loglike_micro]...
    = loglike_compute(data_macro, data_micro,...
      ts_micro, num_smooth_draws, num_interp, num_burnin_periods, M_, oo_, options_)

% Compute likelihood for HANK model

load('approximationParameters.mat', 'nz', 'nMeasure', 'nShare');



%% Macro likelihood and simulation smoother

% Determine variables to smooth
smooth_vars_aux1 = cell(nz,nMeasure,nShare);
smooth_vars_aux2 = cell(nz,nMeasure,nShare);
for i_Measure = 1:nMeasure
    for i_z = 1:nz
        for i_Share = 1:nShare
            smooth_vars_aux1{i_z,i_Measure,i_Share} = sprintf('moment_%d_%d_%d', i_z, i_Measure, i_Share);
            smooth_vars_aux2{i_z,i_Measure,i_Share} = sprintf('measureCoefficient_%d_%d_%d', i_z, i_Measure, i_Share);
        end
    end
end
smooth_vars = char([smooth_vars_aux1(:); smooth_vars_aux2(:)]);

% Run mean smoother and compute macro log likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Micro likelihood per period

T_micro = length(ts_micro);
nobs = dataset_.nobs;

% Fix some warning messages during parallelization
options_new.dataset = [];
options_new.initial_date = [];

% Seeds for simulation smoother
rand_seeds = randi(2^32,1,num_smooth_draws);

loglikes_micro = nan(1,num_smooth_draws);
% disp('Micro likelihood...');
% timer = tic;
% 
% parfor i_draw = 1:num_smooth_draws
%     
%     dataset_fake = struct;
%     dataset_fake.nobs = nobs;
%     
%     % Compute smoothing draw
%     the_smooth_draw = simulation_smoother(smooth_means, smooth_vars, num_burnin_periods, rand_seeds(i_draw), ...
%                                           M_new, oo_new, options_new, dataset_fake, dataset_info, xparam1, estim_params_, bayestopt_);
%                                       
% end
% 
% fprintf('Micro likelihood time: %6.1f sec\n\n', toc(timer));


%% Sum log likelihood

% Micro log likelihood
log_max = max(loglikes_micro);
loglike_micro = log_max + log(mean(exp(loglikes_micro-log_max))); % Formula deals with underflow
if isempty(loglike_micro)
    loglike_micro = 0;
end

loglike = loglike_macro + loglike_micro; % Total log likelihood

end