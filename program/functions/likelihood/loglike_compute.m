function [loglike, loglike_macro, loglike_micro]...
    = loglike_compute(data_macro, ...
                      num_burnin_periods, smooth_vars, num_smooth_draws, ...
                      M_, oo_, options_, ...
                      data_micro, ts_micro, param, ...
                      varargin)

% Compute log likelihood for macro+micro data


%% Macro likelihood and mean smoother

% Run mean smoother and compute macro log likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_, num_smooth_draws>0);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Micro likelihood per period

if isempty(data_micro) || isempty(ts_micro)
    loglike_micro = nan;
    loglike = loglike_macro;
    return;
end

T_micro = length(ts_micro);
nobs = dataset_.nobs;

% Fix some warning messages during parallization
options_new.dataset = [];
options_new.initial_date = [];

% Seeds for simulation smoother
rand_seeds = randi(2^32,1,num_smooth_draws);

% Random shocks for simulation smoother
if isempty(varargin)
    sim_shocks = zeros(0,0,num_smooth_draws);
else
    sim_shocks = varargin{1};
end

loglikes_micro = nan(1,num_smooth_draws);
disp('Micro likelihood...');
timer = tic;

% for i_draw = 1:num_smooth_draws
parfor i_draw = 1:num_smooth_draws
    
    rng(rand_seeds(i_draw), 'twister'); % Set RNG
    
    dataset_fake = struct;
    dataset_fake.nobs = nobs;
    
    % Compute smoothing draw
    the_smooth_draw = simulation_smoother(smooth_means, smooth_vars, num_burnin_periods, sim_shocks(:,:,i_draw), ...
                                          M_new, oo_new, options_new, dataset_fake, dataset_info, xparam1, estim_params_, bayestopt_);
    the_smooth_draw_tab = struct2table(the_smooth_draw); % Transform struct to table
    the_smooth_draw_tab = the_smooth_draw_tab(ts_micro,:); % Only retain relevant time periods for micro data

    the_loglikes_micro_draw = nan(1,T_micro);
    
    try % Once numerical issue in one period, no need to go through the remaining periods, but still run the other smooth draws
        
        for it = 1:T_micro
            
            % Likelihood
            the_likes = likelihood_micro(the_smooth_draw_tab(it,:), permute(data_micro(it,:,:), [2 3 1]), param);
            
            % Log likelihood
            the_loglikes_micro_draw_t = log(the_likes);
            the_loglikes_micro_draw(it) = sum(the_loglikes_micro_draw_t);
            
        end
        
        loglikes_micro(i_draw) = sum(the_loglikes_micro_draw);
    
    catch
        
    end
    
    % Print progress
    if mod(i_draw,ceil(num_smooth_draws/50))==0
        offs = floor(50*i_draw/num_smooth_draws);
        fprintf(['%' num2str(offs+3) 'd%s\n'], round(100*i_draw/num_smooth_draws), '%');
    end
    
end

fprintf('Micro likelihood time: %6.1f sec\n\n', toc(timer));


%% Sum log likelihood

% Micro log likelihood
ix_micro = isfinite(loglikes_micro);
if sum(ix_micro) > 0 % as long as there is a draw survive
    loglikes_micro = loglikes_micro(ix_micro);
    log_max = max(loglikes_micro);
    loglike_micro = log_max + log(mean(exp(loglikes_micro-log_max))); % Formula deals with underflow
else
    loglike_micro = nan;
end
loglike = loglike_macro + loglike_micro; % Total log likelihood


end