function [assets,dist] = dist_irf(maxhorz,shock,M_,oo_,options_,mMoments,mParameters)

    % Impulse response of asset distribution

    global nAssetsFine nEpsilon nMeasure vAssetsGridFine;
    
    assets = vAssetsGridFine;
    dist = nan(nAssetsFine,nEpsilon,maxhorz+1);
    
    % First compute steady state distribution
    dist(:,:,1) = aux_dist(mMoments, mParameters(:,2:end));
    
    % Compute IRFs of all endogenous Dynare variables
    the_irf = computeIRF(maxhorz,shock,M_,oo_,options_);

    % Compute asset distribution at each horizon
    for iHorizon = 1:maxhorz
        
        the_mMoments = nan(nEpsilon,nMeasure);
        the_mMeasureCoefficients = nan(nEpsilon,nMeasure);
        
        % Collect moments and distribution coefficients
        for iEpsilon = 1:nEpsilon
            for iMeasure = 1:nMeasure
                the_mMoments(iEpsilon,iMeasure) = ...
                    the_irf(strcmp(M_.endo_names, ...
                                   ['lag_moment_' num2str(iEpsilon) '_' num2str(iMeasure)]), ...
                            iHorizon);
                the_mMeasureCoefficients(iEpsilon,iMeasure) = ...
                    the_irf(strcmp(M_.endo_names, ...
                                   ['measureCoefficient_' num2str(iEpsilon) '_' num2str(iMeasure)]), ...
                            iHorizon);
            end            
        end
        
        % Compute distribution at this horizon
        dist(:,:,iHorizon+1) = aux_dist(the_mMoments, the_mMeasureCoefficients);
        
    end
    
end


%% Auxiliary functions

function irf_out = computeIRF(maxhorz,shock,M_,oo_,options_)

    % Call internal Dynare impulse response computation
    
    dr = resol(0, M_, options_, oo_); % Decision rule
    irf_out = irf(M_, options_, dr, shock, maxhorz, [], [], 1); % IRFs of all endogenous variables
    irf_out = irf_out+dr.ys; % Add back steady state

end

function dist = aux_dist(mMoments, mMeasureCoefficients)

    % Asset distribution

    global nAssetsFine vAssetsGridFine nEpsilon nMeasure aaBar;

    dist = nan(nAssetsFine,nEpsilon);

    for iEpsilon = 1:nEpsilon

        moment = mMoments(iEpsilon,:);
        measureCoefficient = mMeasureCoefficients(iEpsilon,:);
        moment_aux = moment;
        moment_aux(1) = 0;
        g_log = @(a) measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux');
        normalization = integral(@(a) exp(g_log(a)), aaBar, Inf);

        % Compute density away from borrowing constraint
        dist(:,iEpsilon) = exp(g_log(vAssetsGridFine'))/normalization;
        
    end 

end