function [loglike, loglike_macro, loglike_hh]...
    = loglike_compute_indv_param(data_macro, data_hh,...
      ts_hh, num_smooth_draws, num_burnin_periods, M_, oo_, options_)

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
smooth_vars = char([{'w'; 'r'; 'mHat_1' ; 'mHat_2'}; smooth_vars_aux1(:); smooth_vars_aux2(:)]); %char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones

% Run mean smoother and compute macro likelihood
timer = tic;
[loglike_macro, smooth_means, M_new, oo_new, options_new, dataset_, dataset_info, xparam1, estim_params_, bayestopt_] ...
    = likelihood_smoother(data_macro, smooth_vars, M_, oo_, options_);
fprintf('Macro likelihood/smoother time: %6.1f sec\n\n', toc(timer));


%% Household likelihood per period

T_hh = length(ts_hh);
nobs = dataset_.nobs;
n_lambda = 100;

% Make local versions of global variables so Matlab doesn't complain in the parfor loop
aaBar_local = aaBar;
nMeasure_local = nMeasure;
mmu_local = mmu;
ttau_local = ttau;
mu_l_local = mu_l;

% Seeds for simulation smoother
rand_seeds = randi(2^32,1,num_smooth_draws);

loglikes_hh = nan(1,num_smooth_draws);
disp('Individual likelihood...');
timer = tic;

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
%             n_ix = length(ix);
%             the_likes = nan(1,n_ix);

            % probability
            % prepare the parameters
            mHat = the_smooth_draw.(['mHat_' num2str(eepsilon+1)])(t);
            moment = nan(1,nMeasure_local);
            measureCoefficient = nan(1,nMeasure_local);
            for i_Measure = 1:nMeasure_local
                moment(i_Measure) = the_smooth_draw.(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
                measureCoefficient(i_Measure) = the_smooth_draw.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
            end

            moment_aux = moment;
            moment_aux(1) = 0;
            g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure_local)')-moment_aux'));
            normalization = integral(g, aaBar_local, Inf);
%             g_norm = @(a) g(a)/normalization;
            
            c = the_smooth_draw.w(t)*((1-eepsilon)*mmu_local+eepsilon*(1-ttau_local));
            R = 1+the_smooth_draw.r(t);
            bounds_aux = data_hh(it,ix,2)/(c+R*aaBar_local);
            
            % Continuous part
            % tic
            %             for i_ix = 1:n_ix
            %                 the_likes(i_ix) = (1-mHat)/R*...
            %                                   integral(@(lam) g_norm((data_hh(it,ix(i_ix),2)./lam-c)/R) ...
            %                                                   ./lam ...
            %                                                   .*lognpdf(lam,mu_l_local,sqrt(-2*mu_l_local)), ...
            %                                            0, bounds_aux(i_ix));
            %             end
            
            v_lambda = logninv(linspace(0+1/n_lambda/2,1-1/n_lambda/2,n_lambda)'*logncdf(bounds_aux,mu_l_local,sqrt(-2*mu_l_local)),mu_l_local,sqrt(-2*mu_l_local));
            a_aux = (data_hh(it,ix,2)./v_lambda-c)/R;
            g_aux = measureCoefficient(1)*(a_aux-moment(1));
            for i_Measure = 2:nMeasure_local
                g_aux = g_aux+measureCoefficient(i_Measure)*((a_aux-moment(1)).^i_Measure-moment_aux(i_Measure));
            end
            the_likes = (1-mHat)/R*mean((exp(g_aux)/normalization ...
                ./v_lambda)).*logncdf(bounds_aux,mu_l_local,sqrt(-2*mu_l_local));
            % toc
            
            % Add contribution from point mass
            the_likes = the_likes + mHat/(c+R*aaBar_local)*lognpdf(bounds_aux,mu_l_local,sqrt(-2*mu_l_local));

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

log_max = max(loglikes_hh);
loglike_hh = log_max + log(mean(exp(loglikes_hh-log_max))); % Formula deals with underflow

loglike = loglike_macro + loglike_hh;

end