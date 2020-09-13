clear all;

addpath('./functions/plot');


%% Settings

% Results to plot
model_name = 'hh';      % Either 'hh' or 'firm'
subspec = '_N1000';     % For model_name='hh': format '_N1000'; for model_name='firm': format '_trunc90'
plot_liktypes = 1:2;    % Likelihood types to plot
plot_reps = 1:10;       % Repetitions to include in plot (non-existing repetitions are ignored)

% Plot settings
plot_burnin = 1e3;              % Burn-in
acf_lags = 200;                 % No. of ACF lags
colors_postdens = [zeros(1,3); get(0, 'DefaultAxesColorOrder')]; % Posterior density colors for different liktypes
alpha_postdens = 0.5;
plot_fontsize = 12;             % Plot font size
graph_size = [6 6];             % Graph size for trace plot and ACF
graph_size_postdens = [6 3];    % Graph size for trace plot and ACF

% Parameter names
if strcmp(model_name, 'hh')
    plot_param = {'bbeta', 'ssigmaMeas', 'mu_l'};   % Parameter names to plot
    tex_param = {'\beta','\sigma_e','\mu_\lambda'}; % Tex versions of parameter names (in same order as above)
else
    plot_param = {'rrhoProd', 'ssigmaProd'};
    tex_param = {'\rho_\epsilon','\sigma_\epsilon'};
end

% Folders
results_folder = 'results';                         % Stores results
save_folder = fullfile(results_folder, 'plots');    % Saved figures


%% Load MCMC results

n_liktype = length(plot_liktypes);
n_rep = length(plot_reps);

model_filename = cell(n_rep,n_liktype);
model_mcmc = cell(n_rep,n_liktype);
model_param_names = cell(n_rep,n_liktype);
model_nparam = nan(n_rep,n_liktype);
model_params_truth = cell(n_rep,n_liktype);

% Loop over all relevant results
for i_type = 1:n_liktype
    
    for i_rep = 1:n_rep
        
        % Model file name and results
        model_filename{i_rep,i_type} = sprintf('%s%s%s%d%s%02d',model_name,subspec,'_liktype',plot_liktypes(i_type),'_',plot_reps(i_rep));
        the_file = fullfile(results_folder, model_filename{i_rep,i_type});
        if ~isfile(strcat(the_file, '.mat'))
            continue;
        end
        model_mcmc{i_rep,i_type} = load(the_file); %size ~ 1MB each
        
        % Parameter names and true values
        model_param_names{i_rep,i_type} = model_mcmc{i_rep,i_type}.param_names;
        model_nparam(i_rep,i_type) = length(model_param_names{i_rep,i_type});
        model_params_truth{i_rep,i_type} = nan(1,model_nparam(i_rep,i_type));
        for i_param = 1:model_nparam(i_rep,i_type)
            model_params_truth{i_rep,i_type}(i_param) = model_mcmc{i_rep,i_type}.economicParameters_true.(model_param_names{i_rep,i_type}{i_param});
        end
        
    end
    
end

clearvars the_file;


%% ACF, trace plot, and posterior density

mkdir(save_folder);
nparam_all = length(tex_param);

f_postdens_all = figure;

for i_rep = 1:n_rep
    
    f_postdens = figure;
    the_n_liktype = 0;

    for i_type = 1:n_liktype
        
        if isempty(model_mcmc{i_rep,i_type})
            continue;
        end
        
        the_n_liktype = the_n_liktype+1;
        the_nparam = model_nparam(i_rep,i_type); % Number of parameters for this likelihood type
        
        f_trace = figure;
        f_acf = figure;
        
        for i_param = 1:the_nparam
            
            the_draws = model_mcmc{i_rep,i_type}.post_draws(:,i_param); % Posterior draws for parameter
            the_param = find(strcmp(plot_param,model_param_names{i_rep,i_type}(i_param)),1); % Index of parameter in full parameter vector
            the_tex = tex_param{the_param}; % Tex title for parameter
            
            % Trace plot
            figure(f_trace);
            subplot(the_nparam,1,i_param);
            plot(the_draws,'-k');
            xline(plot_burnin,'k--');
            title(the_tex,'FontSize',plot_fontsize,'FontWeight','bold');
            
            % ACF
            figure(f_acf);
            subplot(the_nparam,1,i_param);
            the_acf = autocorr(the_draws(plot_burnin+1:end),acf_lags);
            plot(0:acf_lags,the_acf,'-k','LineWidth',2);
            yline(0, 'k--');
            title(the_tex,'FontSize',plot_fontsize,'FontWeight','bold');
            
            % Posterior density
            figure(f_postdens);
            subplot(nparam_all,1,the_param);
            [the_f,the_xi] = ksdensity(the_draws(plot_burnin+1:end));
            hold on;
            plot(the_xi,the_f,'LineWidth',2,'Color',colors_postdens(i_type,:));
            hold off;
            xline(model_params_truth{i_rep,i_type}(i_param),'k--');
            title(the_tex,'FontSize',plot_fontsize,'FontWeight','bold');
            
            % Posterior density, all repetitions together
            figure(f_postdens_all);
            subplot(nparam_all,1,the_param);
            hold on;
            patchline(the_xi,the_f,'linestyle','-','edgecolor',colors_postdens(i_type,:),...
                      'linewidth',2,'edgealpha',alpha_postdens);
            hold off;
            xline(model_params_truth{i_rep,i_type}(i_param),'k--');
            title(the_tex,'FontSize',plot_fontsize,'FontWeight','bold');
            
        end
        
        % Save trace plot and ACF
        graph_out(f_trace,fullfile(save_folder,strcat(model_filename{i_rep,i_type},'_trace')),graph_size);
        graph_out(f_acf,fullfile(save_folder,strcat(model_filename{i_rep,i_type},'_acf')),graph_size);
        
    end
    
    % Save posterior densities
    if the_n_liktype>0
        graph_out(f_postdens, ...
                  fullfile(save_folder,...
                           strcat(strrep(model_filename{i_rep,1},sprintf('%s%d','_liktype',plot_liktypes(1)),''), ...
                                  '_postdens') ...
                           ), ...
                  graph_size_postdens);
    else
        close(f_postdens);
    end
    
end

% Save posterior densities, all repetitions together
graph_out(f_postdens_all,fullfile(save_folder,strcat(model_filename{i_rep,1}(1:end-11),'postdens')),graph_size_postdens);

clearvars the_*;


%% Consumption policy functions and distribution IRFs

if strcmp(model_name, 'hh')
    addpath([model_name '_model/plot']);
    consumption_policy_fcn;
    dist_irf;
end