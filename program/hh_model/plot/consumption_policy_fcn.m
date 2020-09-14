% Compute posterior of steady state consumption policy function


%% Calibrate parameters, execute initial Dynare processing

run_calib_dynare;


%% Compute SS consumption policy function

if is_run_comput

    % Consumption policy functions under true DGP
    [vassets_truth,vcons_truth] = aux_cons(M_,oo_,options_);

    % Loop over posterior draws

    comput_ind = find(plot_reps==comput_rep); % Index of repetition used for computations below

    ixs = cell(1,n_liktype);
    assets_draw_all = cell(1,n_liktype);
    cons_draw_all = cell(1,n_liktype);

    for i_type = 1:n_liktype

        ixs{i_type} = plot_burnin+1:thin_draw:model_mcmc{comput_ind,i_type}.mcmc_num_iter;
        assets_draw_all{i_type} = nan(nAssets,length(ixs{i_type}));
        cons_draw_all{i_type} = nan(nAssets,nEpsilon,length(ixs{i_type}));

        for i_draw = 1:length(ixs{i_type})

            fprintf('%d%s%d%s%d%s%d\n', i_draw, '/', length(ixs{i_type}), ': Likelihood ', i_type, ', index = ', ixs{i_type}(i_draw));

            % Update parameters and compute SS policy
            update_param(model_mcmc{comput_ind,i_type}.post_draws(ixs{i_type}(i_draw),:), model_param_names{comput_ind,i_type});
            [assets_draw_all{i_type}(:,i_draw),cons_draw_all{i_type}(:,:,i_draw)] = ...
                aux_cons(M_,oo_,options_);
        end 

    end

    % Save results
    cd('../../');
    save(fullfile(save_folder, 'consumption_policy_fcn.mat'), ...
        'vassets_truth','vcons_truth','assets_draw_all','cons_draw_all','ixs');

else
    
    cd('../../');
    load(fullfile(save_folder, 'consumption_policy_fcn.mat'));

end


%% Plot graphs

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
    plot(vassets_truth,vcons_truth(:,iEpsilon),'k-','linewidth',.75);
    hold off;
    
    xlim(xlim_polfct);
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
graph_out(f_polfct,fullfile(save_folder,'consumption_policy_fcn'),graph_size_polfct);

close all

%% Auxiliary function
function [vassets,vcons] = aux_cons(M_,oo_,options_)

    saveParameters;         % Save parameter values to files
    setDynareParameters;    % Update Dynare parameters in model struct
    compute_steady_state;   % Compute steady state

    expectationCoefficient_mat = nan(nEpsilon,nAssets);
    for i_Asset = 1:nAssets
        for i_Epsilon = 1:nEpsilon
            expectationCoefficient_mat(i_Epsilon,i_Asset) = ...
                eval(['expectationCoefficient_' num2str(i_Epsilon) '_' num2str(i_Asset)]);
        end
    end

    expectationPoly_mat = nan(nAssets,nAssets);
    for i_Power = 1:nAssets
        for i_Asset = 1:nAssets
            expectationPoly_mat(i_Asset,i_Power) = ...
                eval(['expectationPoly_' num2str(i_Asset) '_' num2str(i_Power)]);
        end
    end

    vassets = vAssetsGrid;
    vcons = nan(nAssets,nEpsilon);

    for iEpsilon = 1:nEpsilon
        for iAssets = 1 : nAssets
            s = w*(mmu*(1-(iEpsilon-1))+(1-ttau)*(iEpsilon-1))+(1+r)*vAssetsGrid(iAssets);
            vcons(iAssets,iEpsilon) = min(exp(expectationCoefficient_mat(iEpsilon,:)...
                *expectationPoly_mat(iAssets,:)')^(-1/ssigma),s-aaBar);
        end
    end  

end
