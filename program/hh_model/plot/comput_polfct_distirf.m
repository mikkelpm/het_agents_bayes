% Compute posterior of steady state consumption policy function
% and asset distribution IRF


%% Calibrate parameters, execute initial Dynare processing

run_calib_dynare;


%% Compute consumption policy function and asset distribution IRF

% Index of repetition used for computations below
comput_rep_ind = find(plot_reps==comput_rep);

if is_run_comput
    
    disp('Computing...');
    
    global mat_suff;
    mat_suff = '_comput';
    
    % Maximum impulse response horizon
    maxhorz = max(horzs);

    % Consumption policy functions under true DGP
    timer = tic;
    disp('True parameters');
    update_param(model_params_truth{comput_rep_ind,1}, model_param_names{comput_rep_ind,1});
    [M_new,oo_new,options_new,mMoments,mParameters] = newMod(M_,oo_,options_);
    [assets_truth,cons_truth] = cons_polfct(M_new.steady_vars);
    [assets_fine_truth,distirf_truth] = dist_irf(maxhorz,M_new,oo_new,options_new,mMoments,mParameters);

    % Loop over posterior draws
    ixs = cell(1,n_liktype);
    assets_draw_all = cell(1,n_liktype);
    cons_draw_all = cell(1,n_liktype);
    assets_fine_draw_all = cell(1,n_liktype);
    distirf_draw_all = cell(1,n_liktype);

    for i_type = 1:n_liktype

        ixs{i_type} = plot_burnin+1:thin_draw:model_mcmc{comput_rep_ind,i_type}.mcmc_num_iter; % Indices of posterior draws to use
        
        assets_draw_all{i_type} = nan(nAssets,length(ixs{i_type}));
        cons_draw_all{i_type} = nan(nAssets,nEpsilon,length(ixs{i_type}));
        assets_fine_draw_all{i_type} = nan(nAssetsFine,length(ixs{i_type}));
        distirf_draw_all{i_type} = nan(nAssetsFine,nEpsilon,maxhorz+1,length(ixs{i_type}));

        for i_draw = 1:length(ixs{i_type})

            fprintf('%s%d%s%5d (%3d/%3d)\n', 'Likelihood = ', i_type, ': draw = ', ixs{i_type}(i_draw), i_draw, length(ixs{i_type}));

            % Update parameters and Dynare model structures
            update_param(model_mcmc{comput_rep_ind,i_type}.post_draws(ixs{i_type}(i_draw),:), model_param_names{comput_rep_ind,i_type});
            [M_new,oo_new,options_new,mMoments,mParameters] = newMod(M_,oo_,options_);
            
            % Compute steady state consumption policy function
            [assets_draw_all{i_type}(:,i_draw),cons_draw_all{i_type}(:,:,i_draw)] ...
                = cons_polfct(M_new.steady_vars);
            
            % Compute asset distribution IRF
            [assets_fine_draw_all{i_type}(:,i_draw),distirf_draw_all{i_type}(:,:,:,i_draw)] ...
                = dist_irf(maxhorz,M_new,oo_new,options_new,mMoments,mParameters);
            
        end 

    end
    
    elapsed_time = toc(timer);
    disp('Done with computations. Elapsed time (minutes):');
    disp(elapsed_time/60);
    
    clear M_new oo_new options_new mMoments mParameters;

    % Save results
    cd('../../');
    save(fullfile(save_folder, 'comput_polfct_distirf.mat'), ...
        'ixs','assets_truth','cons_truth','assets_fine_truth','distirf_truth', ...
        'assets_draw_all','cons_draw_all','assets_fine_draw_all','distirf_draw_all', ...
        'elapsed_time');

else
    
    cd('../../');
    load(fullfile(save_folder, 'comput_polfct_distirf.mat'));

end


%% Plot consumption policy function

save_name = strrep(model_filename{comput_rep_ind,1},sprintf('%s%d','_liktype',plot_liktypes(1)),''); % Experiment name without "liktype"
ylim_subplots = nan(nEpsilon,2);
f_polfct = figure;

for iEpsilon = nEpsilon:-1:1
    
    subplot(1,nEpsilon,iEpsilon);
    hold on;
    for i_type = 1:n_liktype
        for i_draw = 1:length(ixs{i_type})
            patchline(assets_draw_all{i_type}(:,i_draw),...
                cons_draw_all{i_type}(:,iEpsilon,i_draw),...
                'linestyle','-','edgecolor',colors_polfct(i_type,:),'linewidth',1,'edgealpha',alpha_polfct);
        end
    end
    plot(assets_truth,cons_truth(:,iEpsilon),'k-','linewidth',.75); % Truth
    hold off;
    
    xlim(xlim_assets);
    ylim_subplots(iEpsilon,:) = get(gca,'ylim');
    xlabel('assets');
    ylabel('consumption');
    title(emp_label{iEpsilon},'FontSize',plot_fontsize,'FontWeight','bold');
    
end

% Enforce same ylim
yylim = [min(ylim_subplots(:)) max(ylim_subplots(:))];
for iEpsilon = 1:nEpsilon
    subplot(1,nEpsilon,iEpsilon);
    ylim(yylim);
end
graph_out(f_polfct,fullfile(save_folder,strcat(save_name,'_conspolfct')),graph_size_polfct);


%% Plot asset distribution IRF

numhorz = length(horzs);
f_distirf = figure;

for iEpsilon = nEpsilon:-1:1
    
    subplot(1,nEpsilon,iEpsilon);
    hold on;
    
    for i_horz = 1:numhorz
        
        the_horz = horzs(i_horz);
        the_y_wedge = wedge_param(1)+wedge_param(2)*i_horz;
        
        for i_type = 1:n_liktype
            for i_draw = 1:length(ixs{i_type})
                patchline(assets_fine_draw_all{i_type}(:,i_draw),...
                    distirf_draw_all{i_type}(:,iEpsilon,the_horz+1,i_draw)+the_y_wedge,...
                    'linestyle','-','edgecolor',colors_distirf(i_type,:),'linewidth',1,'edgealpha',alpha_distirf);
            end
        end
        
        plot(assets_fine_truth,distirf_truth(:,iEpsilon,1)+the_y_wedge,'k-','linewidth',.75); % Truth
        
        text(xlim_assets(1)+0.8*diff(xlim_assets),the_y_wedge+wedge_text,['h = ' num2str(the_horz)],...
            'FontSize',plot_fontsize,'FontWeight','bold'); % Horizon label
        
    end
    
    hold off;
    xlim(xlim_assets);
    ylim(ylim_distirf);
    set(gca,'ytick',[]);
    grid on;
    title(emp_label{iEpsilon},'FontSize',plot_fontsize,'FontWeight','bold');
    
end

graph_out(f_distirf,fullfile(save_folder,strcat(save_name,'_distirf')),graph_size_distirf);


%% Auxiliary function

function [M_,oo_,options_,mMoments,mParameters] = newMod(M_,oo_,options_)

    saveParameters;         % Save parameter values to files
    setDynareParameters;    % Update Dynare parameters in model struct
    compute_steady_state;   % Compute steady state

end

