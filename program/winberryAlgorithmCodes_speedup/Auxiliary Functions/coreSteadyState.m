% Computes market clearing capital stock and associated distribution and 
% decision rules in steady state
%
% Thomas Winberry, July 26th, 2016

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end

% variables in array to pass on
var_array{1}=bbeta;
var_array{2}=ssigma;
var_array{3}=aalpha;
var_array{4}=ddelta;
var_array{5}=aaBar;
var_array{6}=aggEmployment;
var_array{7}=mmu;
var_array{8}=ttau;
var_array{9}=mEpsilonTransition;
var_array{10}=vEpsilonGrid;
var_array{11}=nEpsilon;
var_array{12}=nAssets;
var_array{13}=nState;
var_array{14}=assetsMin;
var_array{15}=assetsMax;
var_array{16}=vAssetsGrid;
var_array{17}=mEpsilonGrid;
var_array{18}=mAssetsGrid;
var_array{19}=vAssetsPoly;
var_array{20}=vAssetsPolySquared;
var_array{21}=mEpsilonPrimeGrid;
var_array{22}=maxIterations;
var_array{23}=tolerance;
var_array{24}=dampening;
var_array{25}=vAssetsPolyFine;
var_array{26}=vAssetsGridFine;
var_array{27}=mEpsilonGridFine;
var_array{28}=mAssetsGridFine;
var_array{29}=nAssetsFine;
var_array{30}=nStateFine;
var_array{31}=vAssetsPolyQuadrature;
var_array{32}=vAssetsGridQuadrature;
var_array{33}=mEpsilonGridQuadrature;
var_array{34}=mAssetsGridQuadrature;
var_array{35}=nAssetsQuadrature;
var_array{36}=vQuadratureWeights;
var_array{37}=vEpsilonInvariant;
var_array{38}=nMeasure;
var_array{39}=splineOpt;
var_array{40}=vAssetsPolyBC;

%----------------------------------------------------------------
% Compute initial guess of market-clearing capital stock using
% histogram approximation of distribution, from Young (2010)
%----------------------------------------------------------------

t0 = tic;
fprintf('Computing initial guess from histogram...\n')

% Solve for market clearing capital stock
f = @(capital) computeMCResidualHistogram(capital,var_array);
options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-2); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
[aggregateCapitalInit,err,exitflag] = fsolve(f,1.01*kRepSS,options);

% Return exitflag if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end	

aggregateCapital = aggregateCapitalInit;
if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%----------------------------------------------------------------
% Compute moments of histogram to use as initial guess for parametric family
%----------------------------------------------------------------

% Compute histogram
[~, mHistogram] = computeMCResidualHistogram(aggregateCapital,var_array);

% Compute moments from histogram
mMomentsHistogram = zeros(nEpsilon,nMeasure);
aGridMoments = zeros(nEpsilon,nAssetsQuadrature,nMeasure); % grid for computing PDF

for iEpsilon = 1 : nEpsilon
	
	% First moment (uncentered)
	mMomentsHistogram(iEpsilon,1) = sum(vAssetsGridFine' .* (mHistogram(iEpsilon,:) ./ ...
		sum(mHistogram(iEpsilon,:))));
	aGridMoments(iEpsilon,:,1) = vAssetsGridQuadrature - mMomentsHistogram(iEpsilon,1);
		
	% Higher order moments (centered)
	for iMoment = 2 : nMeasure
		mMomentsHistogram(iEpsilon,iMoment) = sum(((vAssetsGridFine' - mMomentsHistogram(iEpsilon,1)) .^ iMoment) .* ...
			(mHistogram(iEpsilon,:) ./ sum(mHistogram(iEpsilon,:))));
		aGridMoments(iEpsilon,:,iMoment) = (vAssetsGridQuadrature' - mMomentsHistogram(iEpsilon,1)) .^ ...
			iMoment - mMomentsHistogram(iEpsilon,iMoment);
	end	
	
end

% Mass at borrowing constraint
mHatHistogram = [mHistogram(1,1) / sum(mHistogram(1,:));mHistogram(2,1) / sum(mHistogram(2,:))];

%----------------------------------------------------------------
% Compute market-clearing capital stock from parametric family
%----------------------------------------------------------------

t0 = tic; 
fprintf('Compute steady state from parametric family...\n')

% Solve for market clearing capital stock
f = @(capital) computeMCResidualPolynomials(capital,mMomentsHistogram,aGridMoments,mHatHistogram,var_array);
options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-2);
if abs(f(aggregateCapitalInit)) > 1e-4
	[aggregateCapital,err,exitflag] = fsolve(f,aggregateCapitalInit,options);
end

% Return error if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end	

if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%----------------------------------------------------------------
% Compute other objects from steady state
%----------------------------------------------------------------

[~,mCoefficients,mParameters,mMoments,mHat] = ...
    computeMCResidualPolynomials(aggregateCapital,mMomentsHistogram,aGridMoments,mHatHistogram,var_array);