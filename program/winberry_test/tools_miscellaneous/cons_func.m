clear all;
close all;
clc;

addpath('auxiliary_functions/dynare', 'auxiliary_functions/likelihood', 'auxiliary_functions/sim', 'auxiliary_functions/mcmc');

B_draw = 1e3;
N_draw = 9e3;
thin = 10;

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP ttau mu_l ssigmaMeas;
	

for i_model = 1:2
    if i_model == 1
        load mcmc_long;
    else
        load mcmc_nomicro;
    end
    
    cd('./auxiliary_functions/dynare');
    dynare firstOrderDynamics_polynomials noclearall nopathchange; % Run Dynare once to process model file
    assets_draw = nan(N_draw/thin,nAssets);
    cons_draw = nan(N_draw/thin,nAssets,nEpsilon);
    
    for i_draw = 1:thin:N_draw
        if mod(i_draw,100) == 0
            disp(['Model = ' num2str(i_model) ', i_draw = ' num2str(i_draw)])
        end
        bbeta = post_draws(i_draw,1);
        ssigmaMeas = post_draws(i_draw,2);
        saveParameters;         % Save parameter values to files
        setDynareParameters;    % Update Dynare parameters in model struct
        compute_steady_state;   % Compute steady state once and for all
        
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
        
        assets_draw(i_draw,:) = vAssetsGrid;%typo
        
        for iEpsilon = 1:nEpsilon
            for iAssets = 1 : nAssets
                s = w*(mmu*(1-(iEpsilon-1))+(1-ttau)*(iEpsilon-1))+(1+r)*vAssetsGrid(iAssets);
                cons_draw((i_draw-1)/thin+1,iAssets,iEpsilon) = min((exp(expectationCoefficient_mat(iEpsilon,:)*expectationPoly_mat(iAssets,:)')^(-1/ssigma)),s-aaBar);
            end
        end
    end
    
    save(['cons_model' num2str(i_model) '.mat'],'assets_draw','cons_draw')
    
    for iEpsilon = 1:nEpsilon
        F1 = figure;
        
        hold on
        for i_draw = 1:N_draw/thin
            patchline(assets_draw(i_draw,:)',cons_draw(i_draw,:,iEpsilon)','linestyle','-','edgecolor','b','linewidth',2,'edgealpha',.05);%wrong
        end
        
        % export graph
        graph_tag = ['cons_model' num2str(i_model) '_eps' num2str(iEpsilon)];
        graph_size = [6 6];
        SFont = 16;
        lgd = '';
        graph_out(F1,graph_tag,[],SFont,graph_size,lgd);
        savefig([graph_tag '.fig'])
    end
    cd('../../');
end
