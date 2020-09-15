function [assets,cons] = cons_polfct(ss)

    % Steady state consumption policy function

    global nAssets nEpsilon vAssetsGrid aaBar ssigma mmu ttau;
    
    assets = vAssetsGrid;
    
    polynomials = load_mat('polynomials');

    expectationCoefficient_mat = nan(nEpsilon,nAssets);
    for i_Asset = 1:nAssets
        for i_Epsilon = 1:nEpsilon
            expectationCoefficient_mat(i_Epsilon,i_Asset) = ...
                ss.(['expectationCoefficient_' num2str(i_Epsilon) '_' num2str(i_Asset)]);
        end
    end

    expectationPoly_mat = polynomials.vAssetsPoly;

    cons = nan(nAssets,nEpsilon);

    for iEpsilon = 1:nEpsilon
        for iAssets = 1 : nAssets
            s = ss.w*(mmu*(1-(iEpsilon-1))+(1-ttau)*(iEpsilon-1))+(1+ss.r)*vAssetsGrid(iAssets);
            cons(iAssets,iEpsilon) = min(exp(expectationCoefficient_mat(iEpsilon,:)...
                *expectationPoly_mat(iAssets,:)')^(-1/ssigma),s-aaBar);
        end
    end  

end