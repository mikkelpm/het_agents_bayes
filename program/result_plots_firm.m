%% Settings
B_draw = 1000; % Burn-in
N_draw = 9000; % MCMC draws after burn-in
ix_out = B_draw+(1:N_draw);

model_name = 'firm';
likelihood_type = 1;    % =1: Macro + full-info micro; =2: macro + full-info micro, no truncation; =3: macro + moments micro
v_trunc_quant = [0 0.9]; % Micro sample selection: Lower truncation quantile for labor (steady state distribution)

n_model = length(v_trunc_quant);
model_mcmc = cell(1,n_model);
model_filename = cell(1,n_model);

for i_model = 1:n_model
    model_filename{i_model} = ['results/' model_name '_liktype' num2str(likelihood_type)...
        '_trunc' num2str(v_trunc_quant(i_model)*100) '_01'];
    model_mcmc{i_model} = load(model_filename{i_model});
end

tex_param = {'\rho_\epsilon','\sigma_\epsilon'};
n_param = length(tex_param);

%% MCMC diagnostics
acf_lags = 200;

SFont = 12;
graph_size = [6 6];

for i_model = 1:n_model
    % Trace plot
    F1 = figure;
    for j = 1:n_param
        subplot(n_param,1,j);
        plot(model_mcmc{i_model}.post_draws(1:B_draw+N_draw,j),'-k');
        xline(B_draw,'k--');
        title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
    end
    graph_out(F1,[model_filename{i_model} '_traceplot'],graph_size)
    
    % ACF
    F1 = figure;
    for j = 1:n_param
        subplot(n_param,1,j);
        the_acf = autocorr(model_mcmc{i_model}.post_draws(ix_out,j),acf_lags);
        plot(0:acf_lags,the_acf,'-k','LineWidth',2);
        yline(0, 'k--');
        title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
    end
    graph_out(F1,[model_filename{i_model} '_acf'],graph_size)
end
close all

%% Posterior densities
params_truth = [0.53 0.0364];

co = [zeros(1,3); get(0, 'DefaultAxesColorOrder')];
graph_size = [6 3];

for i_graph = 1:n_model
    F1 = figure;
    for j = 1:n_param
        subplot(1,n_param,j);
        hold on
        for i_model = 1:i_graph
            [the_f,the_xi] = ksdensity(model_mcmc{i_model}.post_draws(ix_out,j));
            plot(the_xi,the_f,'LineWidth',2,'Color',co(i_model,:));
        end
        hold off
        xline(params_truth(j),'k--');
        title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
    end
    graph_out(F1,[model_filename{i_graph} '_postdens'],graph_size)
end
close all

