% based on result_plots_hh.m
clear all;

%% Settings
model_name = 'hh';
n_model = 5;
n_rep = 1;

model_lik = cell(n_rep,n_model);
model_filename = cell(n_rep,n_model);

for i_model = 1:n_model
    for i_rep = 1:n_rep
        model_filename{i_rep,i_model} = ['results/' model_name '_likelihood_N1000' ...
            sprintf('%s%d%s%02d','_liktype',i_model,'_',i_rep)];
        model_lik{i_rep,i_model} = load(model_filename{i_rep,i_model}); %size ~ 1MB each 
    end
end

tex_param = {'\beta','\sigma_e','\mu_\lambda'};
n_param = length(tex_param);

%% Posterior densities
params_truth = [0.96 0.02 -0.25];
param_vals_mult = unique([1 linspace(0.75,1.25,51)]); % Multiples of true parameter to compute
n_lik = length(param_vals_mult);

SFont = 12;
co = [zeros(1,3); get(0, 'DefaultAxesColorOrder')];
graph_size = [10 4];

i_rep = 1;

F1 = figure;
for j = 1:n_param
    subplot(1,n_param,j);
    hold on
    for i_model = 1:n_model
        if i_model == 2 && j == n_param
            continue
        end
        the_lik = model_lik{i_rep,i_model}.lik_all((j-1)*n_lik+(1:n_lik),1);
        plot(params_truth(j)*param_vals_mult,the_lik-max(the_lik),...
            'LineWidth',2,'Color',co(i_model,:));
    end
    hold off
    xline(params_truth(j),'k--');
    title(tex_param{j},'FontSize',SFont,'FontWeight','bold');
    
    if j == 1
        h1 = legend({'FI','macro','3 mom','2 mom','1 mom'},...
            'orientation','horizontal','color','none','FontSize',SFont);
        legend('boxoff')
        
        p = get(h1,'Position');
        p(1) = 0.3;
        p(2) = 0;
        set(h1,'Position',p)
    end
end
graph_out(F1,[model_filename{i_rep,1} '_lik'],graph_size)
savefig([model_filename{i_rep,1} '_lik'])

