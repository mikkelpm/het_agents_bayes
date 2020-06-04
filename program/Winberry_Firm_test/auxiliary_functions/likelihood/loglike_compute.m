function [loglike, loglike_macro, loglike_micro]...
    = loglike_compute(data_macro, data_micro, ts_micro, num_smooth_draws, num_interp, num_burnin_periods, M_, oo_, options_)

% Compute likelihood for firm model

global nMeasure ttheta nnu;


%% Macro likelihood and simulation smoother

% Determine variables to smooth
nMeasure_all = (nMeasure+3)*nMeasure/2;
smooth_vars_aux1 = cell(1,nMeasure_all);
% smooth_vars_aux2 = cell(1,nMeasure_all);
for i_Measure = 1:nMeasure_all
    smooth_vars_aux1{i_Measure} = ['lag_moment_' num2str(i_Measure)];
%     smooth_vars_aux2{i_Measure} = ['measureCoefficient_' num2str(i_Measure)];
end
smooth_vars = [{'logWage'; 'aggregateTFP'}; smooth_vars_aux1(:)]; %; smooth_vars_aux2(:)]);

% Run mean smoother and compute macro log likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_, num_smooth_draws>0);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Micro likelihood per period

if isempty(data_micro)
    loglike_micro = nan;
    loglike = loglike_macro;
    return;
end

T_micro = length(ts_micro);
nobs = dataset_.nobs;

% Make local versions of global variables so Matlab doesn't complain in the parfor loop
% nMeasure_local = nMeasure;
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

        N_micro = length(data_micro(it,:));
%         the_loglikes_micro_draw_t = nan(1,N_micro);
        
        % Prepare smoothed draws
        moment = nan(1,nMeasure_all);
%         measureCoefficient = nan(1,nMeasure_all);
        for i_Measure = 1:nMeasure_all
            moment(i_Measure) = the_smooth_draw.(['lag_moment_' num2str(i_Measure)])(t-1);
%             measureCoefficient(i_Measure) = the_smooth_draw.(['measureCoefficient_' num2str(i_Measure)])(t);
        end
        
%         % Compute normalization constant
%         g = @(prod,logk) exp(g_kernel(prod,logk,moment,measureCoefficient,nMeasure_local));
%         lastwarn('');
%         normalization = integral2(g,-Inf,Inf,-Inf,Inf);
%         warnMsg = lastwarn;
%         if ~isempty(warnMsg)
%             disp(measureCoefficient);
%             error('Improper micro density');
%         end
        
        the_c = log(nnu_local)+the_smooth_draw.aggregateTFP(t)-the_smooth_draw.logWage(t);
        the_jacob = 1-nnu_local;
        the_prod = the_jacob*data_micro(it,:,1)-the_c-ttheta_local*data_micro(it,:,2);
        
        the_mean = moment(1:2);
        the_varcov = [moment(3) moment(4); moment(4) moment(5)];
        the_likes = abs(the_jacob)*mvnpdf([the_prod' data_micro(it,:,2)'], the_mean, the_varcov);
        
%         the_likes = the_jacob*g(the_prod,data_micro(it,:,2))/normalization;
        
%         % Convolution
%         data_aux = (1-nnu_local)*data_micro(it,:)...
%             -nnu_local*(log(nnu_local)-log(the_smooth_draw.wage(t)))....
%             -the_smooth_draw.aggregateTFP(t);
%         the_vals = linspace(min(data_aux),max(data_aux),num_interp); % Compute integral at these grid points for log ouput
%         the_ints = zeros(1,num_interp);
%         for i_in=1:num_interp
%             the_ints(i_in) = integral(@(prod) g(prod,(the_vals(i_in)-prod)/ttheta_local), ...
%                 -Inf, Inf)/(ttheta_local*normalization); 
%         end
%         the_likes = interp1(the_vals,the_ints,data_aux,'pchip'); % Cubic interpolation of integral between grid points
        
        % Log likelihood
        the_loglikes_micro_draw_t = log(the_likes);
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