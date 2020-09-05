clear all;

model_name = 'hh';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));

%% Calibrate parameters and set numerical settings
run([model_name '_model/calibrate']);

cd(['./' model_name '_model/dynare']);
saveParameters;

%% Initial Dynare processing
load('firstOrderDynamics_polynomials_results');
check_matlab_path(false);
dynareroot = dynare_config(); % Add Dynare sub-folders to path

%% Consumption policy functions under true DGP
[vassets0,vcons0] = aux_cons(M_,oo_,options_);

%% Full info vs Macro only	
B_draw = 1000; % Burn-in
N_draw = 9000; % MCMC draws after burn-in
thin_draw = 10;

n_model = 2;
assets_draw_all = nan(nAssets,N_draw/thin_draw,n_model);
cons_draw_all = nan(nAssets,nEpsilon,N_draw/thin_draw,n_model);

for i_model = 1:n_model
    load(['../../results/' model_name '_N1000_model' num2str(i_model) '_01']);
    
    for i_draw = ix_out(1:thin_draw:end)
        ii_draw = (i_draw-B_draw-1)/thin_draw+1; 
            % Index after excluding burn-out and thinning
        if mod(i_draw,100) == 0
            disp(['Model ' num2str(i_model) ', i_draw = ' num2str(i_draw)])
        end
        
        bbeta = post_draws(i_draw,1); % Steady state, so mu_l is irrelevant
        ssigmaMeas = post_draws(i_draw,2);
        [assets_draw_all(:,ii_draw,i_model),cons_draw_all(:,:,ii_draw,i_model)] = ...
            aux_cons(M_,oo_,options_);
    end    
end

%% Save results
cd('../../');
save('results/consumption_policy_fcn.mat','vassets0','vcons0',...
    'assets_draw_all','cons_draw_all')

%% Plot graphs
co = [0.7*ones(1,3); get(0, 'DefaultAxesColorOrder')];
SFont = 12;
graph_size = [6 3];
emp_label = {'Unemployed','Employed'};
ylim_subplots = nan(nEpsilon,2);

F1 = figure;
for iEpsilon = 1:nEpsilon
    subplot(1,nEpsilon,3-iEpsilon)
    hold on  
    for i_model = [2 1]
        for i_draw = 1:N_draw/thin_draw
            patchline(assets_draw_all(:,i_draw,i_model),...
                cons_draw_all(:,iEpsilon,i_draw,i_model),...
                'linestyle','-','edgecolor',co(i_model,:),'linewidth',1,'edgealpha',.05);
        end
    end  
    plot(vassets0,vcons0(:,iEpsilon),'k-','linewidth',.75)
    hold off
    
    xlim([0 10])
    ylim_subplots(iEpsilon,:) = get(gca,'ylim');
    xlabel('asset')
    ylabel('consumption')
    title(emp_label{iEpsilon},'FontSize',SFont,'FontWeight','bold');
end

yylim = [min(ylim_subplots(:)) max(ylim_subplots(:))];
for iEpsilon = 1:nEpsilon
    subplot(1,nEpsilon,3-iEpsilon)
    ylim(yylim);
end
graph_name = 'results/consumption_policy_fcn';
graph_out(F1,graph_name,graph_size);

close all

%% Auxiliary function
function [vassets,vcons] = aux_cons(M_,oo_,options_)
saveParameters;         % Save parameter values to files
setDynareParameters;    % Update Dynare parameters in model struct
compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics

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
