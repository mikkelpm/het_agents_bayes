function [loglike, loglike_macro, loglike_micro]...
    = loglike_compute(data_macro, data_micro,...
      ts_micro, num_smooth_draws, num_burnin_periods, M_, oo_, options_)

% Compute likelihood for HANK model

load('economicParameters.mat', 'bbBar');
load('approximationParameters', 'nz', 'nMeasure', 'nShare');
load('grids.mat','vzInvariant', 'vShareGrid', 'vShareFraction');



%% Macro likelihood and simulation smoother

% Determine variables to smooth
smooth_vars_aux1 = cell(nz,nMeasure,nShare);
smooth_vars_aux2 = cell(nz,nMeasure,nShare);
smooth_vars_aux3 = cell(nz,nShare);
for i_z = 1:nz
    for i_Share = 1:nShare
        for i_Measure = 1:nMeasure
            smooth_vars_aux1{i_z,i_Measure,i_Share} = sprintf('moment_%d_%d_%d', i_z, i_Measure, i_Share);
            smooth_vars_aux2{i_z,i_Measure,i_Share} = sprintf('measureCoefficient_%d_%d_%d', i_z, i_Measure, i_Share);
        end
        smooth_vars_aux3{i_z,i_Share} = sprintf('mHat_%d_%d', i_z, i_Share);
    end
end
smooth_vars = char([smooth_vars_aux1(:); smooth_vars_aux2(:); smooth_vars_aux3(:)]);

% Run mean smoother and compute macro log likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Micro likelihood per period

[T_micro,N_micro,~] = size(data_micro);
nobs = dataset_.nobs;

% Local variables for passing to parallel workers
the_nz = nz;
the_nMeasure = nMeasure;
the_nShare = nShare;
the_vzInvariant = vzInvariant;
the_vShareGrid = vShareGrid;
the_vShareFraction = vShareFraction;
the_bbBar = bbBar;

% Fix some warning messages during parallelization
options_new.dataset = [];
options_new.initial_date = [];

% Seeds for simulation smoother
rand_seeds = randi(2^32,1,num_smooth_draws);

loglikes_micro = nan(1,num_smooth_draws);
disp('Micro likelihood...');
timer = tic;

% for i_draw = 1:num_smooth_draws
parfor i_draw = 1:num_smooth_draws
    
    dataset_fake = struct;
    dataset_fake.nobs = nobs;
    
    % Compute smoothing draw
    the_smooth_draw = simulation_smoother(smooth_means, smooth_vars, num_burnin_periods, rand_seeds(i_draw), ...
                                          M_new, oo_new, options_new, dataset_fake, dataset_info, xparam1, estim_params_, bayestopt_);
    
    the_loglikes_t = nan(1,T_micro);
    
    for i_t = 1:T_micro % For each time period with micro data...
        
        the_t = ts_micro(i_t);
        ix0 = (abs(data_micro(i_t,:,2) - the_bbBar) < 1e-8); % Constrained households in this time period
        
        the_loglikes_t_i_z = nan(N_micro,the_nz);
        
        for i_Share = 1:the_nShare % For each household type...
            
            iSh = (data_micro(i_t,:,1) == the_vShareGrid(i_Share)); % Households with the given profit share
            
            for i_z = 1:the_nz
                the_mHat = the_smooth_draw.(sprintf('mHat_%d_%d',i_z,i_Share))(the_t-1);
                the_moment = nan(1,the_nMeasure);
                the_measureCoefficient = nan(1,the_nMeasure);
                for i_Measure = 1:the_nMeasure    
                    the_moment(i_Measure) = the_smooth_draw.(sprintf('moment_%d_%d_%d',i_z,i_Measure,i_Share))(the_t-1);
                    the_measureCoefficient(i_Measure) = the_smooth_draw.(sprintf('measureCoefficient_%d_%d_%d',i_z,i_Measure,i_Share))(the_t);
                end
                moment_aux = the_moment;
                moment_aux(1) = 0;
                the_g_log = @(b) the_measureCoefficient*((b-the_moment(1)).^((1:the_nMeasure)')-moment_aux');
                the_norm = integral(@(b) exp(the_g_log(b)), the_bbBar, Inf);
                
                the_loglikes_t_i_z(iSh & ix0, i_z) = log(the_mHat); % Log likelihood if constrained
                the_loglikes_t_i_z(iSh & ~ix0, i_z) = log(1-the_mHat) - log(the_norm) ...
                                                    + the_g_log(data_micro(i_t,iSh & ~ix0,2)); % Log likelihood if unconstrained
                
            end
            
            the_loglikes_t_i_z(iSh,:) = the_loglikes_t_i_z(iSh,:) + log(the_vShareFraction(i_Share)); % Add probability of household type
           
        end
        
        the_log_max = max(the_loglikes_t_i_z(:));
        the_loglikes_t = the_log_max + log(exp(the_loglikes_t_i_z - the_log_max)*the_vzInvariant); % Total log likelihood for this time point
    
    end
    
    loglikes_micro(i_draw) = sum(the_loglikes_t); % Total log likelihood given smoothing draw
    
    % Print progress
    if mod(i_draw,ceil(num_smooth_draws/50))==0
        offs = floor(50*i_draw/num_smooth_draws);
        fprintf(['%' num2str(offs+3) 'd%s\n'], round(100*i_draw/num_smooth_draws), '%');
    end
    
end

fprintf('Micro likelihood time: %6.1f sec\n\n', toc(timer));


%% Sum log likelihood

% Micro log likelihood
log_max = max(loglikes_micro);
loglike_micro = log_max + log(mean(exp(loglikes_micro-log_max))); % Formula deals with underflow
if isempty(loglike_micro)
    loglike_micro = 0;
end

loglike = loglike_macro + loglike_micro; % Total log likelihood

end