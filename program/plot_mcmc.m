clear all;

addpath(fullfile('functions', 'plot'));


%% Settings

% Results to plot
model_name = 'hh';      % Either 'hh' or 'firm'
subspec = '_N1000';     % For model_name='hh': format '_N1000'; for model_name='firm': format '_trunc90'
if strcmp(model_name,'hh')
    plot_liktypes = 1:2;% Likelihood types to plot
else
    plot_liktypes = 1;  % Likelihood types to plot
end
plot_reps = 1:10;       % Repetitions to include in plot (non-existing repetitions are ignored)

% Reporting settings
plot_burnin = 1e3;              % Burn-in for MCMC chains
acf_lags = 200;                 % No. of ACF lags

% Parameter names
if strcmp(model_name, 'hh')
    plot_param = {'bbeta', 'ssigmaMeas', 'mu_l'};   % Names of parameters to plot
    tex_param = {'\beta','\sigma_e','\mu_\lambda'}; % Tex versions of parameter names (in same order as above)
else
    plot_param = {'rrhoProd', 'ssigmaProd'};
    tex_param = {'\rho_\epsilon','\sigma_\epsilon'};
end

% Plot layout
colors_postdens = [zeros(1,3); get(0, 'DefaultAxesColorOrder')]; % Posterior density colors for different liktypes
alpha_postdens = 0.5;           % Opacity of density curves when plotted on same figure
plot_fontsize = 12;             % Plot font size
graph_size_diagnostic = [6 6];  % Graph size for trace plot and ACF
graph_size_postdens = [6 3];    % Graph size for posterior density plots

% Extra posterior computations (only for 'hh' model)
is_run_dynare = false;          % Process Dynare model?
is_run_comput = true;           % true: run computation; false: load previous computations from file
comput_rep = 1;                 % Single repetition to use for computations (must be included in "plot_reps")
dynare_model = 'firstOrderDynamics_polynomials'; % Dynare model file for 'hh' model
thin_draw = 10;                 % Compute consumption policy or distribution every X-th draw
emp_label = {'Unemployed','Employed'}; % Labels of employment states
xlim_assets = [0 10];           % Limits of assets axis

% Consumption policy function layout
graph_size_polfct = [6 3];      % Graph size
colors_polfct = [0.7*ones(1,3); get(0, 'DefaultAxesColorOrder')]; % Colors for different liktypes
alpha_polfct = 0.05;            % Opacity

% Asset distribution IRF layout
horzs = [0 2 4 8];              % Impulse response horizons to plot
graph_size_distirf = [6 6];     % Graph size
colors_distirf = colors_postdens; % Colors for different liktypes
alpha_distirf = alpha_polfct;   % Opacity
ylim_distirf = [0 .9];          % y-axis limits
wedge_param = [.8 -.2];         % Shift distributions vertically by wedge_param(1)+wedge_param(2)*horizon, in order to display on single plot
wedge_text = 0.1;               % Shift horizon label vertically by this much relative to above wedge

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
            subplot(1,nparam_all,the_param);
            [the_f,the_xi] = ksdensity(the_draws(plot_burnin+1:end));
            plot(the_xi,the_f,'LineWidth',2,'Color',colors_postdens(i_type,:));
            if isempty(get(get(gca,'Title'),'String'))
                hold on; % hold off will be automatic when the figure is closed
                xline(model_params_truth{i_rep,i_type}(i_param),'k--');
                title(the_tex,'FontSize',plot_fontsize,'FontWeight','bold');
            end
            
            % Posterior density, all repetitions together
            figure(f_postdens_all);
            subplot(1,nparam_all,the_param);
            patchline(the_xi,the_f,'linestyle','-','edgecolor',colors_postdens(i_type,:),...
                      'linewidth',2,'edgealpha',alpha_postdens);
            if isempty(get(get(gca,'Title'),'String'))
                hold on; % hold off will be automatic when the figure is closed
                xline(model_params_truth{i_rep,i_type}(i_param),'k--');
                title(the_tex,'FontSize',plot_fontsize,'FontWeight','bold');
            end
            
        end
        
        % Save trace plot and ACF
        graph_out(f_trace,fullfile(save_folder,strcat(model_filename{i_rep,i_type},'_trace')),graph_size_diagnostic);
        graph_out(f_acf,fullfile(save_folder,strcat(model_filename{i_rep,i_type},'_acf')),graph_size_diagnostic);
        
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


%% Consumption policy function and distribution IRF

if strcmp(model_name, 'hh')
    
    addpath('functions');
    addpath(fullfile('functions', 'likelihood'));
    addpath(fullfile('hh_model', 'dynare'));
    addpath(fullfile('hh_model', 'plot'));
    
    comput_polfct_distirf;
    
end