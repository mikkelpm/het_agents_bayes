clear all;

model_name = 'hh';
n_model = 5;
N_micro = 1e3;          % Number of households per non-missing time period
n_rep = 10*(N_micro==1e3)+1*(N_micro==1e4);

% Load data once for parameters
load(['results/' model_name sprintf('%s%d%s%02d', '_likelihood_N', N_micro, '_', 1)]);

model_lik = cell(n_rep,1);
model_filename = cell(n_rep,1);

for i_rep = 1:n_rep
    model_filename{i_rep} = [model_name sprintf('%s%d%s%02d', '_likelihood_N',...
        N_micro, '_', i_rep)]; % Suffix string for all saved .mat files
    model_lik{i_rep} = load(['results/' model_filename{i_rep}]);
end

tex_param = {'\beta','\sigma_e','\mu_\lambda'};

%% Posterior densities
SFont = 12;
co = [zeros(1,3); get(0, 'DefaultAxesColorOrder')];

v_layout = {'bymodel','byrep'};
n_layout = length(v_layout);
v_xlim = {'symmetric','3mompeak'};
n_xlim = length(v_xlim);
v_ylim = {'orig','minus10'};

aux = [0 cumsum(len_lik)];
ix_lik = repmat([aux(1:end-1)'+1 aux(2:end)'],[1 1 n_xlim]);
for i_param = 1:n_param
    if ~isnan(param_vals_include(i_param,1))
        ix_lik(i_param,1,1) = ix_lik(i_param,1,1)+param_vals_num_coarse;
    end
    if ~isnan(param_vals_include(i_param,2))
        ix_lik(i_param,2,1) = ix_lik(i_param,2,1)-param_vals_num_coarse;
    end
end

for i_xlim = 1:n_xlim  
    all_y = cell(n_rep,n_model,n_param);
    for i_rep = 1:n_model
        for i_model = 1:n_model
            for j = 1:n_param
                ix_aux = ix_lik(j,1,i_xlim):ix_lik(j,2,i_xlim);
                the_lik = model_lik{i_rep}.lik_all(ix_aux,1,i_model);
                all_y{i_rep,i_model,j} = the_lik-max(the_lik);
            end
        end
    end
    
    yylim = nan(n_param,2);
    yylim(:,2) = 0;
    for j = 1:n_param
        aux = cell2mat(all_y(:,:,j));
        yylim(j,1) = min(aux(:));
    end

    for i_layout = 1:n_layout 
        if N_micro == 1e4 && i_layout == 1
            continue
        end
        
        subfolder_tag = ['results/layout' v_layout{i_layout} '_xlim' v_xlim{i_xlim} '_ylim'];
        
        if i_layout == 1
            graph_size = [6 1.7];
            for i_model = 1:n_model
                F1 = figure;
                for j = 1:n_param
                    %                 if ismember(i_model,[2 5]) && j == n_param
                    %                     continue
                    %                 end
                    ix_aux = ix_lik(j,1,i_xlim):ix_lik(j,2,i_xlim);
                    
                    subplot(1,n_param,j);
                    hold on
                    for i_rep = 1:n_rep
                        patchline(lik_grid(ix_aux,j),all_y{i_rep,i_model,j}',...
                            'linestyle','-','edgecolor',co(i_model,:),...
                            'linewidth',2,'edgealpha',.5);
                    end
                    hold off
                    ylim(yylim(j,:))
                    xline(params_truth(j),'k--');
                end
                graph_out(F1,[subfolder_tag v_ylim{1} '/' model_filename{n_rep}...
                    '_model' num2str(i_model)],graph_size)
                
                for j = 1:n_param
                    subplot(1,n_param,j)
                    ylim([-10 0])
                end
                graph_out(F1,[subfolder_tag v_ylim{2} '/' model_filename{n_rep}...
                    '_model' num2str(i_model)],graph_size)
            end
        else
            graph_size = [6 3];
            for i_rep = 1:n_rep
                F1 = figure;
                for j = 1:n_param
                    ix_aux = ix_lik(j,1,i_xlim):ix_lik(j,2,i_xlim);
                    
                    subplot(1,n_param,j);
                    hold on
                    for i_model = 1:n_model
                        %                 if ismember(i_model,[2 5]) && j == n_param
                        %                     continue
                        %                 end
                        plot(lik_grid(ix_aux,j),all_y{i_rep,i_model,j},...
                            'LineWidth',2,'Color',co(i_model,:));
                    end
                    hold off
                    ylim(yylim(j,:))
                    xline(params_truth(j),'k--');
                end
                graph_out(F1,[subfolder_tag v_ylim{1} '/' model_filename{i_rep}],graph_size)
                
                for j = 1:n_param
                    subplot(1,n_param,j)
                    ylim([-10 0])
                end
                graph_out(F1,[subfolder_tag v_ylim{2} '/' model_filename{i_rep}],graph_size)
            end
        end     
    end
end