function mCoefficientsNew = updateCoefficients_polynomials(mCoefficients,var_array)

% Updates polynomial coefficients approximating the conditional expectation function in steady state
% 
% Inputs
%   (1) mCoefficients: nEpsilon x nAssets matrix, storing previous iteration's coefficients
%
% Outputs
%   (1) mCoefficientsNew: nEpsilon x nAssets matrix, storing updated coefficients
% 
% Thomas Winberry, January 19, 2016

% Declare global variables
bbeta=var_array{1};
ssigma=var_array{2};
aalpha=var_array{3};
ddelta=var_array{4};
aaBar=var_array{5};
aggEmployment=var_array{6};
mmu=var_array{7};
ttau=var_array{8};
mEpsilonTransition=var_array{9};
vEpsilonGrid=var_array{10};
nEpsilon=var_array{11};
nAssets=var_array{12};
nState=var_array{13};
assetsMin=var_array{14};
assetsMax=var_array{15};
vAssetsGrid=var_array{16};
mEpsilonGrid=var_array{17};
mAssetsGrid=var_array{18};
vAssetsPoly=var_array{19};
vAssetsPolySquared=var_array{20};
mEpsilonPrimeGrid=var_array{21};
maxIterations=var_array{22};
tolerance=var_array{23};
dampening=var_array{24};
vAssetsPolyFine=var_array{25};
vAssetsGridFine=var_array{26};
mEpsilonGridFine=var_array{27};
mAssetsGridFine=var_array{28};
nAssetsFine=var_array{29};
nStateFine=var_array{30};
vAssetsPolyQuadrature=var_array{31};
vAssetsGridQuadrature=var_array{32};
mEpsilonGridQuadrature=var_array{33};
mAssetsGridQuadrature=var_array{34};
nAssetsQuadrature=var_array{35};
vQuadratureWeights=var_array{36};
vEpsilonInvariant=var_array{37};
nMeasure=var_array{38};
splineOpt=var_array{39};
vAssetsPolyBC=var_array{40};
r=var_array{41};
w=var_array{42};

%---------------------------------------------------------------
% Compute current period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectation = exp(mCoefficients * vAssetsPoly');

% Compute target saving
mAssetsPrimeStar = w * (mmu * (1 - mEpsilonGrid) + (1 - ttau) * mEpsilonGrid) + (1 + r) * mAssetsGrid - ...
	(mConditionalExpectation .^ (-1 / ssigma));

% Compute actual saving
mAssetsPrime = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,nAssets));
mAssetsPrimeGrid = repmat(reshape(mAssetsPrime,1,nState),[nEpsilon 1]);

% Compute next period's polynomials
mAssetsPrimeZeros = scaleDown(mAssetsPrime,assetsMin,assetsMax);
mPolyAssetsPrime = computeChebyshev(nAssets,reshape(mAssetsPrimeZeros,nState,1));

%---------------------------------------------------------------
% Compute next period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectationPrime = exp(mCoefficients * mPolyAssetsPrime');

% Compute target saving
mAssetsPrimePrimeStar = w * (mmu * (1 - mEpsilonPrimeGrid) + (1 - ttau) * mEpsilonPrimeGrid) + (1 + r) * mAssetsPrimeGrid - ...
	(mConditionalExpectationPrime .^ (-1 / ssigma));

% Compute actual savings
mAssetsPrimePrimeGrid = max(mAssetsPrimePrimeStar,aaBar*ones(nEpsilon,nEpsilon*nAssets));

%---------------------------------------------------------------
% Update conditional expectation function
%---------------------------------------------------------------

% Compute new conditional expectation function
mConsumptionPrime = w * (mmu * (1 - mEpsilonPrimeGrid) + (1 - ttau) * mEpsilonPrimeGrid) + (1 + r) * ...
	mAssetsPrimeGrid - mAssetsPrimePrimeGrid;
aConditionalExpectationTilde = reshape(bbeta * mEpsilonTransition * ((1 + r) * (mConsumptionPrime .^ (-ssigma))),...
	nEpsilon,nEpsilon,nAssets);

% Extract the relevant entries
mConditionalExpectation = zeros(nEpsilon,nAssets);
for iEpsilon = 1:nEpsilon
	mConditionalExpectation(iEpsilon,:) = aConditionalExpectationTilde(iEpsilon,iEpsilon,:);
end

% Update the coefficients
mCoefficientsNew = zeros(nEpsilon,nAssets);
for iEpsilon = 1:nEpsilon
	vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * log(mConditionalExpectation(iEpsilon,:))),2);
	mCoefficientsNew(iEpsilon,:) = (vCoefficients ./ vAssetsPolySquared)';
end

