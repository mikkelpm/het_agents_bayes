clear all;

addpath(fullfile('functions', 'plot'));


%% Settings

% Results to plot
model_name = 'hh';      
subspec = '_N1000';     % Format '_N1000'
plot_liktypes = 1:5;    % Likelihood types to plot
plot_reps = 1:10;       % Repetitions to include in plot (non-existing repetitions are ignored)

% Parameter names
plot_param = {'bbeta', 'ssigmaMeas', 'mu_l'};   % Names of parameters to plot
tex_param = {'\beta','\sigma_e','\mu_\lambda'}; % Tex versions of parameter names (in same order as above)

% Plot layout
colors_lik = [zeros(1,3); 1 0 0; get(0, 'DefaultAxesColorOrder')]; % Posterior density colors for different liktypes
colors_lik(4,:) = [];           % Delete color(s) too similar
alpha_lik = 0.5;                % Opacity
plot_fontsize = 12;             % Plot font size
graph_size_rep = [6 2.3];       % Graph size for plots by repetitions
graph_size_liktype = [6 1.7];   % Graph size for plots by likelihood types
ylim_lik = [-10 0];             % Limits of likelihood axis
legend_lik = {'FI','macro','3 mom','2 mom','1 mom'}; % Likelihood type names in legend
shrink_rep = 0.05;               % Shrink the plots by repetitions due to legend

% Folders
results_folder = 'results';                         % Stores results
save_folder = fullfile(results_folder, 'plots');    % Saved figures

% For each parameter, which likelihood type(s) cannot identify the parameter?
nonid_param_lik = {[],[],[2 5]};

%% Load likelihood results

n_rep = length(plot_reps);

model_filename = cell(n_rep,1);
model_lik = cell(n_rep,1);

% Loop over all relevant results
for i_rep = 1:n_rep
    
    % Model file name and results
    model_filename{i_rep} = sprintf('%s%s%s%s%02d',model_name,'_likelihood',...
                                    subspec,'_',plot_reps(i_rep));
    the_file = fullfile(results_folder, model_filename{i_rep});
    if ~isfile(strcat(the_file, '.mat'))
        continue;
    end
    model_lik{i_rep} = load(the_file); 
    
end

clearvars the_file;

% Additional variables passed from simulation
var_list = {'len_lik','lik_grid','params_truth'};
for i_var = 1:length(var_list)
    eval(sprintf('%s = model_lik{1}.%s;',var_list{i_var},var_list{i_var}));
end


%% Likelihood plots

n_liktype = length(plot_liktypes);
n_param = length(tex_param);

all_param = cell(n_param,1);             % Vector for each parameter
all_lik = cell(n_rep,n_liktype,n_param); % Standardized likelhoods (max = 0)
xlim_param = nan(n_param,2,n_rep);       % Limits of parameter axes

aux = [0 cumsum(len_lik)];
aux_lim = [aux(1:end-1)'+1 aux(2:end)'];

for i_param = 1:n_param
    
    ix_aux = aux_lim(i_param,1):aux_lim(i_param,2);
    all_param{i_param} = lik_grid(ix_aux,i_param);
    
    for i_rep = 1:n_rep
        
        xlim_aux = nan(n_liktype,2);
        
        for i_type = 1:n_liktype 
        
            if ismember(plot_liktypes(i_type),nonid_param_lik{i_param})
                continue
            end
            
            the_lik = model_lik{i_rep}.lik_all(ix_aux,1,i_type);
            all_lik{i_rep,i_type,i_param} = the_lik-max(the_lik);  
            
            aux = all_param{i_param}(all_lik{i_rep,i_type,i_param}>=ylim_lik(1));
            xlim_aux(i_type,:) = [min(aux) max(aux)];
            
        end
        
        xlim_param(i_param,:,i_rep) = [min(xlim_aux(:,1)) max(xlim_aux(:,2))];
        
    end 
    
end

% Lsyout 1: by repetitions
for i_rep = 1:n_rep
    
    f_rep = figure;
    
    for i_param = 1:n_param
        
        subplot(1,n_param,i_param); 
        hold on
        for i_type = 1:n_liktype
            
            if ismember(plot_liktypes(i_type),nonid_param_lik{i_param})
                continue
            end
            
            plot(all_param{i_param},all_lik{i_rep,i_type,i_param},...
                'LineWidth',1,'Color',colors_lik(i_type,:));
            
        end
        hold off
        
        xlim(xlim_param(i_param,:,i_rep))
        ylim(ylim_lik)
        xline(params_truth(i_param),'k--');
        title(tex_param{i_param},'FontSize',plot_fontsize,'FontWeight','bold');
        
        subplot_position = get(gca,'Position');
        subplot_position(2) = subplot_position(2)+shrink_rep;
        subplot_position(4) = subplot_position(4)-shrink_rep;
        set(gca,'Position',subplot_position)
        
        if i_param == 1
            
            legend_handle = legend(legend_lik(plot_liktypes),... 
                'orientation','horizontal','color','none');
            legend('boxoff')
            legend_handle.Position(1) = 0.5 - legend_handle.Position(3)/2;
            legend_handle.Position(2) = 0;
            
        end
        
    end
    
    graph_out(f_rep,fullfile(save_folder,model_filename{i_rep}),graph_size_rep);
        
end

% Layout 2: by likelihood types
xlim_param_all = [min(xlim_param(:,1,:),[],3) max(xlim_param(:,2,:),[],3)];

for i_type = 1:n_liktype
    
    f_type = figure;
    
    for i_param = 1:n_param
        
        if ismember(plot_liktypes(i_type),nonid_param_lik{i_param})
            continue
        end
        
        subplot(1,n_param,i_param);
        hold on
        for i_rep = 1:n_rep
            patchline(all_param{i_param},all_lik{i_rep,i_type,i_param},...
                'linestyle','-','edgecolor',colors_lik(i_type,:),...
                'linewidth',1,'edgealpha',alpha_lik);
        end
        hold off
        
        xlim(xlim_param_all(i_param,:))
        ylim(ylim_lik)
        xline(params_truth(i_param),'k--');

    end

    graph_out(f_type,fullfile(save_folder,sprintf('%s%s%d',...
        model_filename{1}(1:end-2),'liktype',plot_liktypes(i_type))),...
        graph_size_liktype);
    
end
