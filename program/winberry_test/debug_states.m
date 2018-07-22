clear all
close all
clc

addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim');
debug_folder = 'auxiliary_functions\debug_netinc_1e4\';
load([debug_folder 'loglike_space_96.mat']);

n_smooth_vars = size(smooth_vars,1);

e_aux = data_hh(:,:,1);
e_aux = e_aux(:);
a_aux = data_hh(:,:,3);
a_aux = a_aux(:);
n_aplot = 100;
v_aplot = [linspace(0,max(a_aux(e_aux==0)),n_aplot);
    linspace(0,max(a_aux(e_aux==1)),n_aplot)];

%% smooth draws
all_states = nan(T_hh, n_smooth_vars, num_smooth_draws);
all_fplot = nan(n_aplot,num_smooth_draws,2,T_hh);

for i_draw = 1:num_smooth_draws
    % parfor i_draw = 1:num_smooth_draws
    
    dataset_fake = struct;
    dataset_fake.nobs = nobs;
    
    % Compute smoothing draw
    the_smooth_draw = simulation_smoother(smooth_means, smooth_vars, num_burnin_periods, rand_seeds(i_draw), ...
        M_new, oo_new, options_new, dataset_fake, dataset_info, xparam1, estim_params_, bayestopt_);
    for iter_j=1:n_smooth_vars
        the_var = deblank(smooth_vars(iter_j,:));
        all_states(:,iter_j,i_draw) = the_smooth_draw.(the_var)(ts_hh);
    end
    
    
    for it = 1:T_hh
        
        t = ts_hh(it);
        
        for eepsilon = 0:1
            
            % probability
            % prepare the parameters
            moment = nan(1,nMeasure_local);
            measureCoefficient = nan(1,nMeasure_local);
            for i_Measure = 1:nMeasure_local
                moment(i_Measure) = the_smooth_draw.(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
                measureCoefficient(i_Measure) = the_smooth_draw.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
            end
            
            % Compute normalization constant
            moment_aux = moment;
            moment_aux(1) = 0;
            g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure_local)')-moment_aux'));
            lastwarn('');
            normalization = integral(g, aaBar_local, Inf);
            warnMsg = lastwarn;
            if ~isempty(warnMsg)
                disp(measureCoefficient);
                error('Improper asset density');
            end
            
            all_fplot(:,i_draw,eepsilon+1,it) = g(v_aplot(eepsilon+1,:))/normalization;
            
        end
    end
end

%% true states
benchmark = load([debug_folder 'simul.mat']);
benchmark_states = nan(T_hh, n_smooth_vars);
benchmark_fplot = nan(n_aplot,2,T_hh);

for iter_j=1:n_smooth_vars
    the_var = deblank(smooth_vars(iter_j,:));
    benchmark_states(:,iter_j) = benchmark.(the_var)(ts_hh);
end

dataset_fake = struct;
dataset_fake.nobs = nobs;
the_loglikes_hh_draw = nan(1,T_hh);

for it = 1:T_hh
    
    t = ts_hh(it);
    
    the_loglikes_hh_draw_t = nan(1,length(data_hh(it,:,1)));
    
    for eepsilon = 0:1
        
        ix = find(data_hh(it,:,1)==eepsilon);
        
        % probability
        % prepare the parameters
        mHat = benchmark.(['mHat_' num2str(eepsilon+1)])(t);            
        moment = nan(1,nMeasure_local);
        measureCoefficient = nan(1,nMeasure_local);
        for i_Measure = 1:nMeasure_local
            moment(i_Measure) = benchmark.(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
            measureCoefficient(i_Measure) = benchmark.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
        end
        
        % Compute normalization constant
        moment_aux = moment;
        moment_aux(1) = 0;
        g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure_local)')-moment_aux'));
        lastwarn('');
        normalization = integral(g, aaBar_local, Inf);
        warnMsg = lastwarn;
        if ~isempty(warnMsg)
            disp(measureCoefficient);
            error('Improper asset density');
        end
        
        % Continuous part
        c = benchmark.w(t)*((1-eepsilon)*mmu_local+eepsilon*(1-ttau_local));
        R = benchmark.r(t);
        if R<=0
            warning('%s%8.4f', 'R=', R);
        end
        bounds_aux = data_hh(it,ix,2)/(c+R*aaBar_local);
        
        v_lambda = logninv(linspace(0+1/n_lambda/2,1-1/n_lambda/2,n_lambda)'*logncdf(bounds_aux,mu_l_local,sqrt(-2*mu_l_local)),mu_l_local,sqrt(-2*mu_l_local));
        a_aux = (data_hh(it,ix,2)./v_lambda-c)/R;
        g_aux = measureCoefficient(1)*(a_aux-moment(1));
        for i_Measure = 2:nMeasure_local
            g_aux = g_aux+measureCoefficient(i_Measure)*((a_aux-moment(1)).^i_Measure-moment_aux(i_Measure));
        end
        the_likes = (1-mHat)/R*mean((exp(g_aux)/normalization ...
            ./v_lambda)).*logncdf(bounds_aux,mu_l_local,sqrt(-2*mu_l_local));
        
        % Point mass
        the_likes = max(the_likes + mHat/(c+R*aaBar_local)*lognpdf(bounds_aux,mu_l_local,sqrt(-2*mu_l_local)),eps);
        
        the_loglikes_hh_draw_t(ix) = log(the_likes);
        
        benchmark_fplot(:,eepsilon+1,it) = g(v_aplot(eepsilon+1,:))/normalization;
            
    end
    the_loglikes_hh_draw(it) = sum(the_loglikes_hh_draw_t);
end

benchmark_ll = sum(the_loglikes_hh_draw)

% %% plot state
% h1 = figure;
% set(h1,'Units','Inches','position',[0 0 10 4])
% set(h1,'PaperOrientation','portrait');
% set(h1,'PaperPositionMode','Auto');
% set(h1,'Clipping','off')
% SFont = 16;
% 
% for iter_j=1:n_smooth_vars
%     the_var = deblank(smooth_vars(iter_j,:));
%     for it = 1:T_hh
%         subplot(1,T_hh,it)
%         [ff, xx] = ksdensity(squeeze(all_states(it,iter_j,:)));
%         plot(xx,ff,'b-','linewidth',2)
%         yylim = [0 max(ff)*1.1];
%         hold on
%         plot(benchmark_states(it,iter_j)*ones(1,2),yylim,'r-','linewidth',2);
%         hold off
%         xlim([min(xx(1),benchmark_states(it,iter_j)) max(xx(end),benchmark_states(it,iter_j))])
%         ylim(yylim)
%         title(['t=' num2str(ts_hh(it))])
%         set(gca,'FontSize',SFont,'FontWeight','bold')
%     end
%     print(h1,[debug_folder 'state_' the_var],'-dpng');
% end
% 
%% plot dist

for it = 1:T_hh
    for eepsilon = 0:1
        h1 = figure;
        set(h1,'Units','Inches','position',[0 0 6 6])
        set(h1,'PaperOrientation','portrait');
        set(h1,'PaperPositionMode','Auto');
        set(h1,'Clipping','off')
        SFont = 16;
        hold on
        for i_draw = 1:num_smooth_draws
            patchline(v_aplot(eepsilon+1,:),all_fplot(:,i_draw,eepsilon+1,it),'linestyle','-','edgecolor',[.7 .7 .7],'linewidth',1,'edgealpha',0.05);
        end
        plot(v_aplot(eepsilon+1,:),benchmark_fplot(:,eepsilon+1,it),'k-','linewidth',2);
        [ff,xx] = ksdensity(data_hh(it,data_hh(it,:,1)==eepsilon,3));
        plot(xx,ff,'b-','linewidth',2)
        hold off
        xlim(v_aplot(eepsilon+1,[1 end]))
        set(gca,'FontSize',SFont,'FontWeight','bold')
        print(h1,[debug_folder 't' num2str(ts_hh(it)) 'e' num2str(eepsilon)],'-dpng');

    end
end
