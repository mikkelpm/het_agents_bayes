%% Settings
B_draw = 1000; % Burn-in
N_draw = 9000; % MCMC draws after burn-in
ix_out = B_draw+(1:N_draw);

% Alternative setup for n repetitions
N_draw_alt = 4000; % MCMC draws after burn-in
ix_out_alt = B_draw+(1:N_draw_alt);

model_name = 'hh';
n_model = 5;
n_rep = 10;

model_mcmc = cell(n_rep,n_model);
model_filename = cell(n_rep,n_model);

for i_model = 1:n_model
    for i_rep = 1:n_rep
        model_filename{i_rep,i_model} = ['results/' model_name '_N1000' ...
            sprintf('%s%d%s%02d','_liktype',i_model,'_',i_rep)];
        model_mcmc{i_rep,i_model} = load(model_filename{i_rep,i_model}); %size ~ 1MB each 
    end
end

% =========================================================================
% N = 100
model_filename100 = ['results/' model_name '_N100_model1_01'];
model_mcmc100 = load(model_filename100);

tex_param = {'\beta','\sigma_e','\mu_\lambda'};
n_param = length(tex_param);

%% MCMC diagnostics
acf_lags = 200;

SFont = 12;
graph_size = [6 6];

for i_model = 1:n_model
    for i_rep = 1:n_rep
        if i_model <= 2 && i_rep == 1
            N_draw_aux = N_draw;
        else
            N_draw_aux = N_draw_alt;
        end
            
        % Trace plot
        F1 = figure;
        for j = 1:n_param
            if i_model == 2 && j == n_param
                continue
            end
            
            subplot(n_param,1,j);
            plot(model_mcmc{i_rep,i_model}.post_draws(1:B_draw+N_draw_aux,j),'-k');
            xline(B_draw,'k--');
            title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
        end
        graph_out(F1,[model_filename{i_rep,i_model} '_traceplot'],graph_size)
        
        % ACF
        if i_rep == 1
            F1 = figure;
            for j = 1:n_param
                if i_model == 2 && j == n_param
                    continue
                end
                
                subplot(n_param,1,j);
                the_acf = autocorr(model_mcmc{i_rep,i_model}.post_draws(ix_out_aux,j),acf_lags);
                plot(0:acf_lags,the_acf,'-k','LineWidth',2);
                yline(0, 'k--');
                title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
            end
            graph_out(F1,[model_filename{i_rep,i_model} '_acf'],graph_size)
        end
    end
end

% =========================================================================
% N = 100
% Trace plot
F1 = figure;
for j = 1:n_param
    subplot(n_param,1,j);
    plot(model_mcmc100.post_draws(1:B_draw+N_draw,j),'-k');
    xline(B_draw,'k--');
    title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
end
graph_out(F1,[model_filename100 '_traceplot'],graph_size)

% ACF
F1 = figure;
for j = 1:n_param
    subplot(n_param,1,j);
    the_acf = autocorr(model_mcmc100.post_draws(ix_out,j),acf_lags);
    plot(0:acf_lags,the_acf,'-k','LineWidth',2);
    yline(0, 'k--');
    title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
end
graph_out(F1,[model_filename100 '_acf'],graph_size)

close all

%% Posterior densities
params_truth = [0.96 0.02 -0.25];

co = [zeros(1,3); get(0, 'DefaultAxesColorOrder')];
graph_size = [6 3];

% Full info vs Macro only
i_rep = 1;

F1 = figure;
for j = 1:n_param
    subplot(1,n_param,j);
    hold on
    for i_model = 1:2
        if i_model == 2 && j == n_param
            continue
        end
        [the_f,the_xi] = ksdensity(model_mcmc{i_rep,i_model}.post_draws(ix_out,j));
        plot(the_xi,the_f,'LineWidth',2,'Color',co(i_model,:));
    end
    hold off
    xline(params_truth(j),'k--');
    title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
end
graph_out(F1,[model_filename{i_rep,1} '_postdens'],graph_size)

% =========================================================================
% N = 100
F1 = figure;
for j = 1:n_param
    subplot(1,n_param,j);
    hold on
    for i_model = 1:2
        if i_model == 2 && j == n_param
            continue
        end
        
        if i_model == 1
            [the_f,the_xi] = ksdensity(model_mcmc00.post_draws(ix_out,j));
        else
            [the_f,the_xi] = ksdensity(model_mcmc{i_rep,i_model}.post_draws(ix_out,j));
        end
        plot(the_xi,the_f,'LineWidth',2,'Color',co(i_model,:));
    end
    hold off
    xline(params_truth(j),'k--');
    title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
end
graph_out(F1,[model_filename00 '_postdens'],graph_size)

% =========================================================================
% All 5 models x 10 repetitions
ix_rep = 1:n_rep;
all_xi = cell(n_rep,n_model,n_param);
all_f = cell(n_rep,n_model,n_param);
for i_model = 1:n_model  
    for j = 1:n_param
        if i_model == 2 && j == n_param
            continue
        end
        for i_rep = ix_rep
            [all_f{i_rep,i_model,j},all_xi{i_rep,i_model,j}] = ...
                ksdensity(model_mcmc{i_rep,i_model}.post_draws(ix_out_alt,j));
        end
    end
end

xxlim = nan(n_param,2);
yylim = nan(n_param,2); 
yylim(:,1) = 0;
for j = 1:n_param
    aux = cell2mat(all_xi(:,:,j));
    xxlim(j,:) = [min(aux(:)) max(aux(:))];
    aux = cell2mat(all_f(:,:,j));
    yylim(j,2) = max(aux(:));
end

for i_model = 1:n_model  
    F1 = figure;
    if i_model == 1
        graph_size = [6 1.8]; % due to subtitle
        position_inch = nan(n_param,4); % inner position in inches,
                                        % so all graphs are compatible across models
    else
        graph_size = [6 1.7];
    end
    
    for j = 1:n_param
        if i_model == 2 && j == n_param
            continue
        end
        
        subplot(1,n_param,j);
        hold on
        for i_rep = ix_rep
            patchline(all_xi{i_rep,i_model,j}',all_f{i_rep,i_model,j}','linestyle','-',...
                'edgecolor',co(i_model,:),'linewidth',2,'edgealpha',.5);
        end
        hold off
        
        xlim(xxlim(j,:))
        ylim(yylim(j,:))
        xline(params_truth(j),'k--');
        
        if i_model == 1
            title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
            position_inch(j,:) = get(gca,'Position').*repmat(graph_size,2);
        else
            set(gca,'Units','Inches','Position',position_inch(j,:))
        end
    end
    
    graph_out(F1,[model_filename{n_rep,i_model} '_postdens'],graph_size)
end

close all

%% Consumption policy functions and distribution IRFs
addpath([model_name '_model/plot']);
consumption_policy_fcn;
dist_irf;