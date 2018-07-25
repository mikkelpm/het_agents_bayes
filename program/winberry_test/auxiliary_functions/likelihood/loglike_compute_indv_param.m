function [loglike, loglike_macro, loglike_hh]...
    = loglike_compute_indv_param(data_macro, data_hh,...
      ts_hh, num_smooth_draws, num_interp, num_burnin_periods, M_, oo_, options_)

% Compute likelihood for Krusell-Smith model

global aaBar nEpsilon nMeasure mmu ttau mu_l;


%% Macro likelihood and simulation smoother

% Determine variables to smooth
smooth_vars_aux1 = cell(nEpsilon,nMeasure);
smooth_vars_aux2 = cell(nEpsilon,nMeasure);
for i_Measure = 1:nMeasure
    for i_Epsilon = 1:nEpsilon
        smooth_vars_aux1{i_Epsilon,i_Measure} = ['moment_' num2str(i_Epsilon) '_' num2str(i_Measure)];
        smooth_vars_aux2{i_Epsilon,i_Measure} = ['measureCoefficient_' num2str(i_Epsilon) '_' num2str(i_Measure)];
    end
end
smooth_vars = char([{'w'; 'r'; 'mHat_1' ; 'mHat_2'}; smooth_vars_aux1(:); smooth_vars_aux2(:)]);

% Run mean smoother and compute macro log likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Household likelihood per period

T_hh = length(ts_hh);
nobs = dataset_.nobs;

% Make local versions of global variables so Matlab doesn't complain in the parfor loop
aaBar_local = aaBar;
nMeasure_local = nMeasure;
mmu_local = mmu;
ttau_local = ttau;
mu_l_local = mu_l;

% Fix some warning messages during parallization
options_new.dataset = [];
options_new.initial_date = [];

% Seeds for simulation smoother
rand_seeds = randi(2^32,1,num_smooth_draws);

loglikes_hh = nan(1,num_smooth_draws);
disp('Individual likelihood...');
timer = tic;

% for i_draw = 1:num_smooth_draws
parfor i_draw = 1:num_smooth_draws
    
    dataset_fake = struct;
    dataset_fake.nobs = nobs;
    
    % Compute smoothing draw
    the_smooth_draw = simulation_smoother(smooth_means, smooth_vars, num_burnin_periods, rand_seeds(i_draw), ...
                                          M_new, oo_new, options_new, dataset_fake, dataset_info, xparam1, estim_params_, bayestopt_);
    
    the_loglikes_hh_draw = nan(1,T_hh);

    for it = 1:T_hh
    
        t = ts_hh(it);

        the_loglikes_hh_draw_t = nan(1,length(data_hh(it,:,1)));

        for eepsilon = 0:1

            ix = find(data_hh(it,:,1)==eepsilon);

            % Prepare smoothed draws
            mHat = the_smooth_draw.(['mHat_' num2str(eepsilon+1)])(t);
            moment = nan(1,nMeasure_local);
            measureCoefficient = nan(1,nMeasure_local);
            for i_Measure = 1:nMeasure_local
                moment(i_Measure) = the_smooth_draw.(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
                measureCoefficient(i_Measure) = the_smooth_draw.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
            end
            
            % Compute normalization constant
            moment_aux = moment;
            moment_aux(1) = 0;
            g_log = @(a) measureCoefficient*((a-moment(1)).^((1:nMeasure_local)')-moment_aux');
            lastwarn('');
            normalization = integral(@(a) exp(g_log(a)), aaBar_local, Inf);
            warnMsg = lastwarn;
            if ~isempty(warnMsg)
                disp(measureCoefficient);
                error('Improper asset density');
            end
            
            % Continuous part
            c = the_smooth_draw.w(t)*((1-eepsilon)*mmu_local+eepsilon*(1-ttau_local));
            R = the_smooth_draw.r(t);
            if R<=0
                warning('%s%8.4f', 'R=', R);
            end
            
            the_sigma2 = -2*mu_l_local;
            
            the_vals = linspace(min(log(data_hh(it,ix,2))),max(log(data_hh(it,ix,2))),num_interp); % Compute integral at these grid points for log income
            the_ints = zeros(1,num_interp);
            for i_in=1:num_interp
                the_ints(i_in) = integral(@(a) exp(g_log(a) ...
                                                   -(0.5/the_sigma2)*((the_vals(i_in)-mu_l_local)-log(c+R*a)).^2 ...
                                                   ), ...
                                          aaBar_local, Inf); % Compute integral at given grid point
            end
            the_likes = interp1(the_vals,the_ints,log(data_hh(it,ix,2)),'pchip'); % Cubic interpolation of integral between grid points
            the_likes = (the_likes./data_hh(it,ix,2)) * ((1-mHat)/(normalization*sqrt(2*pi*the_sigma2))); % Likelihood for continuous part
            
            % Point mass part
            the_likes = max(the_likes ...
                            + (mHat/sqrt(2*pi*the_sigma2))*exp(-0.5/the_sigma2*(log(data_hh(it,ix,2))-log(c+R*aaBar_local)-mu_l_local).^2)./data_hh(it,ix,2), ...
                            eps);
            
            % Log likelihood
            the_loglikes_hh_draw_t(ix) = log(the_likes);

        end
        
        the_loglikes_hh_draw(it) = sum(the_loglikes_hh_draw_t);

    end
    
    loglikes_hh(i_draw) = sum(the_loglikes_hh_draw);
    
    % Print progress
    if mod(i_draw,ceil(num_smooth_draws/50))==0
        offs = floor(50*i_draw/num_smooth_draws);
        fprintf(['%' num2str(offs+3) 'd%s\n'], round(100*i_draw/num_smooth_draws), '%');
    end
    
end

fprintf('Individual likelihood time: %6.1f sec\n\n', toc(timer));


%% Sum log likelihood

% Micro log likelihood
log_max = max(loglikes_hh);
loglike_hh = log_max + log(mean(exp(loglikes_hh-log_max))); % Formula deals with underflow
if isempty(loglike_hh)
    loglike_hh = 0;
end

loglike = loglike_macro + loglike_hh; % Total log likelihood

end