function [loglike, loglike_macro, loglike_micro]...
    = loglike_compute(data_macro, data_micro,...
      ts_micro, num_smooth_draws, num_interp, num_burnin_periods, M_, oo_, options_)

% Compute likelihood for Krusell-Smith model

global nMeasure ttheta nnu;


%% Macro likelihood and simulation smoother

% Determine variables to smooth
nMeasure_all = (nMeasure+3)*nMeasure/2;
smooth_vars_aux1 = cell(1,nMeasure_all);
smooth_vars_aux2 = cell(1,nMeasure_all);
for i_Measure = 1:nMeasure_all
    smooth_vars_aux1{i_Measure} = ['moment_' num2str(i_Measure)];
    smooth_vars_aux2{i_Measure} = ['measureCoefficient_' num2str(i_Measure)];
end
smooth_vars = char([{'wage'; 'aggregateTFP'}; smooth_vars_aux1(:); smooth_vars_aux2(:)]);

% Run mean smoother and compute macro log likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Micro likelihood per period

T_micro = length(ts_micro);
nobs = dataset_.nobs;

% Make local versions of global variables so Matlab doesn't complain in the parfor loop
nMeasure_local = nMeasure;
ttheta_local = ttheta;
nnu_local = nnu;

% Fix some warning messages during parallization
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
    
    the_loglikes_micro_draw = nan(1,T_micro);

    for it = 1:T_micro
    
        t = ts_micro(it);

        n_micro = length(data_micro(it,:));
        the_loglikes_micro_draw_t = nan(1,n_micro);
        
        % Prepare smoothed draws
        mHat = the_smooth_draw.(['mHat_' num2str(eepsilon+1)])(t);
        moment = nan(1,nMeasure_all);
        measureCoefficient = nan(1,nMeasure_all);
        for i_Measure = 1:nMeasure_all
            moment(i_Measure) = the_smooth_draw.(['moment_' num2str(i_Measure)])(t-1);
            measureCoefficient(i_Measure) = the_smooth_draw.(['measureCoefficient_' num2str(i_Measure)])(t-1);
        end
        
        % Compute normalization constant
        g = @(prod,logk) exp(g_kernel(prod,logk,moment,measureCoefficient));
        lastwarn('');
        normalization = integral2(g,-Inf,Inf,-Inf,Inf);
        warnMsg = lastwarn;
        if ~isempty(warnMsg)
            disp(measureCoefficient);
            error('Improper asset density');
        end
        
        % Convolution
        data_aux = (1-nnu_local)*data_micro(it,:)...
            -nnu_local*(log(nnu_local)-log(the_smooth_draw.wage(t)))....
            -the_smooth_draw.aggregateTFP(t);
        the_likes = zeros(1,n_micro);
        for i_micro=1:n_micro
            the_likes(i_micro) = integral(@(prod) g(prod,(data_aux(i_micro)-prod)/ttheta_local)/normalization, ...
                -inf, Inf); 
        end
        
        % Log likelihood
        the_loglikes_micro_draw_t(ix) = log(the_likes);
        the_loglikes_micro_draw(it) = sum(the_loglikes_micro_draw_t);

    end
    
    loglikes_micro(i_draw) = sum(the_loglikes_micro_draw);
    
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