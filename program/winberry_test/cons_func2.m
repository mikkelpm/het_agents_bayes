global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
    mmu rrhoTFP ssigmaTFP ttau mu_l ssigmaMeas;

cd('./auxiliary_functions/dynare');

bbeta = .96;
ssigmaMeas = .02;
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

vassets = vAssetsGrid;
vcons= nan(nAssets,nEpsilon);

for iEpsilon = 1:nEpsilon
    for iAssets = 1 : nAssets
        s = w*(mmu*(1-(iEpsilon-1))+(1-ttau)*(iEpsilon-1))+(1+r)*vAssetsGrid(iAssets);
        vcons(iAssets,iEpsilon) = min((exp(expectationCoefficient_mat(iEpsilon,:)*expectationPoly_mat(iAssets,:)')^(-1/ssigma)),s-aaBar);
    end
end


for i_model = 1:2
    load(['cons_model' num2str(i_model) '.mat'])
    
    if i_model == 1
        assets_draw_all = assets_draw(1001:10:end,:);
        cons_draw_all = cons_draw(101:end,:,:);
    else
        assets_draw_all(:,:,2) = assets_draw(1001:10:end,:);
        cons_draw_all(:,:,:,2) = cons_draw(101:end,:,:);
    end
end

cd('../../');
cons_draw_same = nan(3,25,2,2);
for i4 = 1:2
    for i3 = 1:2
        cons_draw_aux = nan(800,25);
        for i1 = 1:800
            cons_draw_aux(i1,:) = interp1(assets_draw_all(i1,:,i3),cons_draw_all(i1,:,i3,i4),vassets);
        end
        cons_draw_same(:,:,i3,i4) = quantile(cons_draw_aux,[.005 .5 .995]);
    end
end

yylim = [min(cons_draw_same(:))-0.05 max(cons_draw_same(:))+0.05];
for iEpsilon = 1:2
    F1 = figure;
    
    hold on
    plot(vassets,cons_draw_same(2,:,iEpsilon,1),'k-','linewidth',2)
%     plot(vassets,cons_draw_same([1 3],:,iEpsilon,1),'k--','linewidth',1)
    jbfill(vassets',cons_draw_same(3,:,iEpsilon,1),cons_draw_same(1,:,iEpsilon,1),'k','k',[],0.3)
    plot(vassets,cons_draw_same(2,:,iEpsilon,2),'r-','linewidth',2)
%     plot(vassets,cons_draw_same([1 3],:,iEpsilon,2),'r--','linewidth',1)
    jbfill(vassets',cons_draw_same(3,:,iEpsilon,2),cons_draw_same(1,:,iEpsilon,2),'r','r',[],0.3)
    plot(vassets,vcons(:,iEpsilon),'b-','linewidth',1)
    hold off
    
    xlim(vassets([1 20]))
    ylim([min(reshape(cons_draw_same(:,1:20,iEpsilon,:),[],1))-0.01 max(reshape(cons_draw_same(:,1:20,iEpsilon,:),[],1))+0.01])
%     ylim(yylim)
    
    % export graph
    graph_tag = ['cons_eps' num2str(iEpsilon)];
    graph_size = [4 4];
    SFont = 16;
    lgd = '';
    graph_out(F1,graph_tag,[],SFont,graph_size,lgd);
    savefig([graph_tag '.fig'])
end

for iEpsilon = 1:2
    F1 = figure;
    
    hold on
    
    for i_draw = 1:800
        patchline(assets_draw_all(i_draw,:,2)',cons_draw_all(i_draw,:,iEpsilon,2)','linestyle','-','edgecolor','r','linewidth',1,'edgealpha',.05);
    end
    
    for i_draw = 1:800
        patchline(assets_draw_all(i_draw,:,1)',cons_draw_all(i_draw,:,iEpsilon,1)','linestyle','-','edgecolor',[.7 .7 .7],'linewidth',1,'edgealpha',.05);
    end
%     plot(vassets,cons_draw_same(2,:,iEpsilon,1),'k-','linewidth',2)
%     plot(vassets,cons_draw_same([1 3],:,iEpsilon,1),'k--','linewidth',1)
%     jbfill(vassets',cons_draw_same(3,:,iEpsilon,1),cons_draw_same(1,:,iEpsilon,1),'k','k',[],0.3)
%     plot(vassets,cons_draw_same(2,:,iEpsilon,2),'r-','linewidth',2)
%     plot(vassets,cons_draw_same([1 3],:,iEpsilon,2),'r--','linewidth',1)
%     jbfill(vassets',cons_draw_same(3,:,iEpsilon,2),cons_draw_same(1,:,iEpsilon,2),'r','r',[],0.3)
    plot(vassets,vcons(:,iEpsilon),'b-','linewidth',.75)
    hold off
    
%     xlim(vassets([1 25]))
    xlim([vassets(1) 6.5])
%     ylim([min(reshape(cons_draw_same(:,1:20,iEpsilon,:),[],1))-0.01 max(reshape(cons_draw_same(:,1:20,iEpsilon,:),[],1))+0.01])
    xlabel('asset')
    ylabel('consumption')
%     ylim(yylim)
    
    % export graph
    graph_tag = ['cons_eps' num2str(iEpsilon) '_v2'];
    graph_size = [4 4];
    SFont = 14;
    lgd = '';
    graph_out(F1,graph_tag,[],SFont,graph_size,lgd);
    savefig([graph_tag '.fig'])
end





