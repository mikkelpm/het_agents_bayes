%% Settings
B_draw = 1000; % Burn-in
N_draw = 9000; % MCMC draws after burn-in
ix_out = B_draw+(1:N_draw);

% Alternative setup for n repetitions
N_draw_alt = 4000; % MCMC draws after burn-in
ix_out_alt = B_draw+(1:N_draw_alt);

model_name = 'firm';
likelihood_type = 1;    % =1: Macro + full-info micro; =2: macro + full-info micro, no truncation; =3: macro + moments micro
v_trunc_quant = [0 0.9]; % Micro sample selection: Lower truncation quantile for labor (steady state distribution)
n_model = length(v_trunc_quant);
n_rep = 10;

model_mcmc = cell(n_rep,n_model);
model_filename = cell(n_rep,n_model);

for i_model = 1:n_model
    for i_rep = 1:n_rep
        model_filename{i_rep,i_model} = ['results/' model_name ...
            '_liktype' num2str(likelihood_type)...
            '_trunc' num2str(v_trunc_quant(i_model)*100) ...
            '_' sprintf('%02d',i_rep)];
        try
            model_mcmc{i_rep,i_model} = load(model_filename{i_rep,i_model});
        catch
        end
    end
end

tex_param = {'\rho_\epsilon','\sigma_\epsilon'};
n_param = length(tex_param);

%% MCMC diagnostics
acf_lags = 200;

SFont = 12;
graph_size = [6 6];

for i_model = 1:n_model
    if i_model == 1
        v_rep = 1:n_rep;
    else
        v_rep = 1;
    end
    
    for i_rep = v_rep
        if i_rep == 1
            N_draw_aux = N_draw;
        else
            N_draw_aux = N_draw_alt;
        end
        
        % Trace plot
        F1 = figure;
        for j = 1:n_param
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
                subplot(n_param,1,j);
                the_acf = autocorr(model_mcmc{i_rep,i_model}.post_draws(ix_out,j),acf_lags);
                plot(0:acf_lags,the_acf,'-k','LineWidth',2);
                yline(0, 'k--');
                title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
            end
            graph_out(F1,[model_filename{i_rep,i_model} '_acf'],graph_size)
        end
    end
end
close all

%% Posterior densities 
params_truth = [0.53 0.0364];

co = [zeros(1,3); get(0, 'DefaultAxesColorOrder')];
graph_size = [6 3];

% Repetition 1
i_rep = 1;

for i_graph = 1:n_model
    F1 = figure;
    for j = 1:n_param
        subplot(1,n_param,j);
        hold on
        for i_model = 1:i_graph
            [the_f,the_xi] = ksdensity(model_mcmc{i_rep,i_model}.post_draws(ix_out,j));
            plot(the_xi,the_f,'LineWidth',2,'Color',co(i_model,:));
        end
        hold off
        xline(params_truth(j),'k--');
        title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
    end
    graph_out(F1,[model_filename{i_rep,i_graph} '_postdens'],graph_size)
end

% All 10 repetitions, no truncations
i_model = 1;
ix_rep = 1:n_rep;

F1 = figure;
for j = 1:n_param
    subplot(1,n_param,j);
    hold on
    for i_rep = ix_rep
        [the_f,the_xi] = ksdensity(model_mcmc{i_rep,i_model}.post_draws(ix_out_alt,j));
        patchline(the_xi,the_f,'linestyle','-','edgecolor',co(i_model,:),...
            'linewidth',2,'edgealpha',.5);
    end
    hold off
    xline(params_truth(j),'k--');
    
    title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
end
graph_out(F1,[model_filename{n_rep,i_model} '_postdens'],graph_size)

close all

