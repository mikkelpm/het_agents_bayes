clear all;

model_name = 'hh';

addpath(genpath('./functions'));
addpath(genpath(['./' model_name '_model/auxiliary_functions']));

is_run_dynare = true;   % Process Dynare model?

%% Calibrate parameters and set numerical settings
run([model_name '_model/calibrate']);

cd(['./' model_name '_model/dynare']);
saveParameters;

%% Initial Dynare processing
if is_run_dynare
    dynare firstOrderDynamics_polynomials_distirf noclearall nopathchange;
        % Run Dynare once to process model file
else
    load('firstOrderDynamics_polynomials_distirf_results');
    check_matlab_path(false);
    dynareroot = dynare_config(); % Add Dynare sub-folders to path
end

%% Steady state distributions under true DGP
[vassets0,vdist0] = aux_distss(M_,oo_,options_);

%% Full info vs Macro only	
B_draw = 1000; % Burn-in
N_draw = 9000; % MCMC draws after burn-in
thin_draw = 10;

n_model = 2;
nHorizon = 8;
assets_draw_all = nan(nAssets,N_draw/thin_draw,n_model);
distirf_draw_all = nan(nAssets,nEpsilon,nHorizon+1,N_draw/thin_draw,n_model);
v_discard = cell(1,n_model);

for i_model = 1:n_model
    load(['../../results/' model_name '_N1000_model' num2str(i_model) '_01']);
    v_discard{i_model} = [];
    
    for i_draw = ix_out(1:thin_draw:end)
        ii_draw = (i_draw-B_draw-1)/thin_draw+1; 
            % Index after excluding burn-out and thinning
        if mod(i_draw,100) == 0
            disp(['Model ' num2str(i_model) ', i_draw = ' num2str(i_draw)])
        end
        
        bbeta = post_draws(i_draw,1); 
        ssigmaMeas = post_draws(i_draw,2);
        mu_l = post_draws(i_draw,3);
        try
            [~,distirf_draw_all(:,:,1,ii_draw,i_model)] = aux_distss(M_,oo_,options_);
            [assets_draw_all(:,ii_draw,i_model),...
                distirf_draw_all(:,:,2:end,ii_draw,i_model)] = aux_distirf(M_,oo_,options_);
        catch
            v_discard{i_model} = [v_discard{i_model} i_draw];
            continue
        end
    end    
end

%% Save results
cd('../../');
save('results/dist_irf.mat','vassets0','vdist0','assets_draw_all','distirf_draw_all')

%% Plot graphs
co = [0.7*ones(1,3); get(0, 'DefaultAxesColorOrder')];
SFont = 12;
graph_size = [6 6];
emp_label = {'Unemployed','Employed'};
vHorizon_out = [0 2 4 8];
nHorizon_out = length(vHorizon_out);

F1 = figure;
for iEpsilon = 1:nEpsilon
    subplot(1,nEpsilon,3-iEpsilon) 
    
    hold on
    for iHorizon_out = 1:nHorizon_out
        iHorizon = vHorizon_out(iHorizon_out);
        y_wedge = .8-.2*iHorizon_out;
        
        for i_model = [2 1]
            for i_draw = 1:N_draw/thin_draw
                patchline(assets_draw_all(:,i_draw,i_model),...
                    distirf_draw_all(:,iEpsilon,iHorizon+1,i_draw,i_model)+y_wedge,...
                    'linestyle','-','edgecolor',co(i_model,:),'linewidth',1,'edgealpha',.05);
            end
        end
        plot(vassets0,vdist0(:,iEpsilon)+y_wedge,'k-','linewidth',.75)        
    end
    hold off
    
    xlim([0 10])
    ylim([0 .9])
    set(gca,'ytick',[])
    grid on    
    
    for iHorizon_out = 1:nHorizon_out
        iHorizon = vHorizon_out(iHorizon_out);
        text(8,.9-.2*iHorizon_out,['h = ' num2str(iHorizon)],...
            'FontSize',SFont,'FontWeight','bold')
    end
    title(emp_label{iEpsilon},'FontSize',SFont,'FontWeight','bold');
end

graph_name = 'results/dist_irf';
graph_out(F1,graph_name,graph_size);

close all

%% Auxiliary functions
% Distributions at steady states
function [vassets,vdist] = aux_distss(M_,oo_,options_)
saveParameters;         % Save parameter values to files
setDynareParameters;    % Update Dynare parameters in model struct
compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics

vassets = vAssetsGridFine;
vdist = nan(nAssetsFine,nEpsilon);

for iEpsilon = 1:nEpsilon
    moment = mMoments(iEpsilon,:);
    measureCoefficient = mParameters(iEpsilon,2:nMeasure+1);
    moment_aux = moment;
    moment_aux(1) = 0;
    g_log = @(a) measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux');
    normalization = integral(@(a) exp(g_log(a)), aaBar, Inf);
    
    % Compute density away from borrowing constraint
    vdist(:,iEpsilon) = exp(g_log(vAssetsGridFine'))/normalization;
end 
end

% Distribution IRFs
function [vassets,vdist] = aux_distirf(M_,oo_,options_)
saveParameters;         % Save parameter values to files
firstOrderDynamics_polynomials;

vassets = vAssetsGridFine;
vdist = nan(nAssetsFine,nEpsilon,nHorizon);

for iEpsilon = 1:nEpsilon
    for iHorizon = 1:nHorizon
        moment = nan(1,nMeasure);
        for iMeasure = 1:nMeasure
            moment(1,iMeasure) = ...
                eval(['oo_.irfs.lag_moment_' num2str(iEpsilon) '_' num2str(iMeasure)...
                '_aggregateTFPShock(iHorizon)']);
        end
        measureCoefficient = nan(1,nMeasure);
        for iMeasure = 1:nMeasure
            measureCoefficient(1,iMeasure) = ...
                eval(['oo_.irfs.measureCoefficient_' num2str(iEpsilon) '_' num2str(iMeasure)...
                '_aggregateTFPShock(iHorizon)']);
        end
        
        moment_aux = moment;
        moment_aux(1) = 0;
        g_log = @(a) measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux');
        normalization = integral(@(a) exp(g_log(a)), aaBar, Inf);
        
        % Compute density away from borrowing constraint
        vdist(:,iEpsilon,iHorizon) = exp(g_log(vAssetsGridFine'))/normalization;        
    end
end
end
