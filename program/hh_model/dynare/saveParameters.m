% Generated grids/polynomials and save parameter values


%% Call parameters (the next set of commands will overwrite some)

setParameters;


%% Save parameters in .mat files to import into Dynare 

if splineOpt == 0	% if using polynomials to approximate individual decisions

	% Economic parameters
	save_mat economicParameters bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP ssigmaMeas mu_l mEpsilonTransition
		
	% Approximation parameters
	save_mat approximationParameters nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save_mat grids vAssetsGridZeros vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
		
	% Polynomials
	save_mat polynomials vAssetsPoly vAssetsPolySquared vAssetsPolyFine vAssetsPolyQuadrature vAssetsPolyBC
	
else	% if using splines to approximate individual decisions

	% Economic parameters
	save_mat economicParameters bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save_mat approximationParameters nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save_mat grids vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
	
end


