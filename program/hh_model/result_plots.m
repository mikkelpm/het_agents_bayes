%% Settings

B_draw = 1000;
N_draw = 8000;
ix_out = B_draw+(1:N_draw);
acf_lags = 200;

n_model = 4;
model_mcmc = cell(1,n_model);
model_mcmc{1} = load('temp/mcmc_long.mat','post_draws');
model_mcmc{2} = load('temp/mcmc_long2.mat','post_draws');
% model_mcmc{3} = load('temp/mcmc_long3.mat','post_draws');
model_mcmc{4} = load('temp/mcmc_nomicro.mat','post_draws');

l_model = {'FIMicro', '3rdMoment', '2ndMoment', 'Macro'};
tex_param = {'\beta','\sigma_e','\mu_\lambda'};
n_param = length(tex_param);

ix_model = [1 2 4];

%% MCMC diagnostics

for i_model = ix_model
    % Trace plot
    F1 = figure('Units', 'normalize', 'Position', [0.2 0.1 0.6 0.8]);
    for j=1:n_param
        if i_model == 4 && j == n_param
            continue
        end
        subplot(n_param,1,j);
        plot(model_mcmc{i_model}.post_draws(1:B_draw+N_draw,j), '-k');
        xline(B_draw,'k--');
        title(tex_param{j}, 'FontSize', 14);
    end
    savefig(F1,['temp/' l_model{i_model} '_traceplot'])
    
    % ACF
    F1 = figure('Units', 'normalize', 'Position', [0.2 0.1 0.6 0.8]);
    for j=1:n_param
        if i_model == 4 && j == n_param
            continue
        end
        subplot(n_param,1,j);
        the_acf = autocorr(model_mcmc{i_model}.post_draws(ix_out,j), acf_lags);
        plot(0:acf_lags, the_acf, '-k', 'LineWidth', 2);
        yline(0, 'k--');
        title(tex_param{j}, 'FontSize', 14);
    end
    savefig(F1,['temp/' l_model{i_model} '_acf'])
end

%% Posterior densities

params_truth = [0.9600    0.0200   -0.2500];

co =[   0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560];

F1 = figure('Units', 'normalize', 'Position', [0.1 0.2 0.8 0.6]);
for j=1:n_param
    subplot(1,n_param,j);
    hold on
    for i_model = ix_model
        if i_model == 4 && j == n_param
            continue
        end
        [the_f,the_xi] = ksdensity(model_mcmc{i_model}.post_draws(ix_out,j));
        plot(the_xi, the_f, 'LineWidth', 2, 'Color', co(i_model,:));
    end
    hold off
    xline(params_truth(j),'k--');
    title(tex_param{j}, 'FontSize', 14);
    
    if j == 1
        h1 = legend(l_model(ix_model),'orientation','horizontal','color','none','FontSize',14);
        legend('boxoff')
        
        p = get(h1,'Position');
        p(1) = 0.3;
        p(2) = 0;
        set(h1,'Position',p)
    end
end

savefig(F1,'temp/post_dens')

