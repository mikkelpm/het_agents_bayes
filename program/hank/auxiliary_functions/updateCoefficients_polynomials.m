function [mCoefficientsNew,mConditionalExpectation] = updateCoefficients_polynomials(mCoefficients,...
    bbeta,ppsi,nnu,bbBar,eepsilon,ttau,vvarthetaT,mzTransition,nz,nAssets,nState,assetsMin,assetsMax,...
    mzGrid,mAssetsGrid,vAssetsPoly,vAssetsPolySquared,mzPrimeGrid,...
    A_SS,w_SS,rr,NN)

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

%---------------------------------------------------------------
% Compute current period's savings policy function
%---------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectation = exp(mCoefficients * vAssetsPoly');

% Compute target saving
chidT_SS = (1/eepsilon + vvarthetaT)*A_SS*NN;
mAssetsPrimeStar = ((1-ttau)*w_SS*mzGrid).^(1+1/nnu).*(mConditionalExpectation/ppsi).^(1/nnu) ...
                    +(1+rr)*mAssetsGrid+chidT_SS-1./mConditionalExpectation;

% Compute actual saving
mAssetsPrime = max(mAssetsPrimeStar,bbBar);
mAssetsPrimeGrid = repmat(reshape(mAssetsPrime,1,nState),[nz 1]);

% Compute next period's polynomials
mAssetsPrimeZeros = min(max(2 * ((mAssetsPrime - assetsMin) / (assetsMax - assetsMin)) - 1,-1),1); 
mPolyAssetsPrime = ones(nState,nAssets);
mPolyAssetsPrime(:,2) = mAssetsPrimeZeros(:);
for iPower = 3:nAssets
	mPolyAssetsPrime(:,iPower) = 2 * mAssetsPrimeZeros(:) .* mPolyAssetsPrime(:,iPower-1) - mPolyAssetsPrime(:,iPower-2);
end

%---------------------------------------------------------------
% Update conditional expectation function
%---------------------------------------------------------------

% Compute new conditional expectation function
% ASSUMES nnu=1!
    
% Compute conditional expectation
mConditionalExpectationPrime = exp(mCoefficients * mPolyAssetsPrime');

% Compute target saving
mAssetsPrimePrimeStar = ((1-ttau)*w_SS*mzPrimeGrid).^(1+1/nnu).*(mConditionalExpectationPrime/ppsi).^(1/nnu) ...
                +(1+rr)*mAssetsPrimeGrid+chidT_SS-1./mConditionalExpectationPrime;

mConstr = (mAssetsPrimePrimeStar < bbBar);

mLaborPrime = (1-ttau)*w_SS*mzPrimeGrid.*mConditionalExpectationPrime/ppsi;
aux = -bbBar + (1+rr)*mAssetsPrimeGrid(mConstr) + chidT_SS;
mLaborPrime(mConstr) = (-aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*mzPrimeGrid(mConstr)).^2/ppsi)) ...
                       ./ (2*(1-ttau)*w_SS*mzPrimeGrid(mConstr));
mConsumptionPrime = (1-ttau)*w_SS*mzPrimeGrid./(ppsi*mLaborPrime);

aConditionalExpectationTilde = reshape(bbeta * mzTransition * ((1 + rr) ./ mConsumptionPrime),...
	nz,nz,nAssets);

% Extract the relevant entries
mConditionalExpectation = zeros(nz,nAssets);
for iz = 1:nz
	mConditionalExpectation(iz,:) = aConditionalExpectationTilde(iz,iz,:);
end

% Update the coefficients
mCoefficientsNew = zeros(nz,nAssets);
for iz = 1:nz
    vCoefficients = vAssetsPoly' * log(mConditionalExpectation(iz,:))';
	mCoefficientsNew(iz,:) = (vCoefficients ./ vAssetsPolySquared)';
end

