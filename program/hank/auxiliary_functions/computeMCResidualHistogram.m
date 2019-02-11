function [residual,mHistogramOptional,mAssetsPrimeOptional,mConsumptionOptional,mLaborOptional,mConditionalExpectationOptional] = ...
    computeMCResidualHistogram(rr,NN,var_array)

% Computes residual of market-clearing condition, using histogram approximation of distribution
% as in Young (2010); used to compute initial guess for parametric family
% 
% Inputs
%   (1) capital: candidate aggregate capital stock
%
% Outputs
%   (1) residual: residual of market clearing condition (to use in root finder)
%   (2) (optional) mHistogramOptional: mHistogram - histogram distribution
%   (3) (optional) mAssetsPrimeOptional: mAssetsPrime
%   (4) (optional) mConsumptionOptional: mConsumption
% 
% Thomas Winberry, July 26th, 2016

% Declare global variables
bbeta = var_array{1};
ppsi = var_array{2};
nnu = var_array{3};
bbBar = var_array{4};
eepsilon = var_array{5};
ttau = var_array{6};
vvarthetaB = var_array{7};
vvarthetaT = var_array{8};
mzTransition = var_array{9};
vzInvariant = var_array{10};
vShareGrid = var_array{11};
vShareFraction = var_array{12};
nz = var_array{13};
nAssets = var_array{14};
nState = var_array{15};
nStateFine = var_array{16};
nAssetsFine = var_array{17};
nAssetsQuadrature = var_array{18};
nMeasure = var_array{19};
nShare = var_array{20};
assetsMin = var_array{21};
assetsMax = var_array{22};
vzGrid = var_array{23};
mzGrid = var_array{24};
mzGridFine = var_array{25};
mzPrimeGrid = var_array{26};
mzGridQuadrature = var_array{27};
vAssetsGridFine = var_array{28};
vAssetsGridQuadrature = var_array{29};
vAssetsPoly = var_array{30};
vAssetsPolySquared = var_array{31};
vAssetsPolyFine = var_array{32};
vAssetsPolyQuadrature = var_array{33};
vAssetsPolyBC = var_array{34};
mAssetsGrid = var_array{35};
mAssetsGridFine = var_array{36};
mAssetsGridQuadrature = var_array{37};
vQuadratureWeights = var_array{38};
maxIterations = var_array{39};
tolerance = var_array{40};
dampening = var_array{41};
numNewton = var_array{42};
A_SS = var_array{43};
w_SS = var_array{44};


mHistogramOptional = nan(nShare,nz,nAssetsFine);
mAssetsPrimeOptional = nan(nShare,nz,nAssetsFine);
mConsumptionOptional = nan(nShare,nz,nAssetsFine);
mLaborOptional = nan(nShare,nz,nAssetsFine);
mConditionalExpectationOptional = nan(nShare,nz,nAssets);
bond_demand = nan(nShare,1);
labor_supply = nan(nShare,1);

for iShare = 1:nShare

    %----------------------------------------------------------------
    % Compute individual decisions
    %----------------------------------------------------------------

    chidT_SS = (1/eepsilon*vShareGrid(iShare) + vvarthetaT)*A_SS*NN; %chi*d+T

    % Initialize coefficients using rule of thumb savings rule that sets b'=b and z'=z
    % ASSUMES nnu=1!
    aux = rr*mAssetsGrid + chidT_SS;
    mGridInit = log(bbeta * (1+rr) * 2 ./ (aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*mzGrid).^2/ppsi)));

    mCoefficients = zeros(nz,nAssets);
    for iz = 1:nz	% interpolate
        vCoefficients = vAssetsPoly' * mGridInit(iz,:)';
        mCoefficients(iz,:) = (vCoefficients ./ vAssetsPolySquared)';
    end

    % Iterate
    err = Inf; iteration = 1;
    while err > tolerance && iteration <= maxIterations

        [mCoefficientsNew,mConditionalExpectation_hist] = updateCoefficients_polynomials(mCoefficients,...
            bbeta,ppsi,nnu,bbBar,eepsilon,ttau,vvarthetaT,mzTransition,nz,nAssets,nState,assetsMin,assetsMax,...
            mzGrid,mAssetsGrid,vAssetsPoly,vAssetsPolySquared,mzPrimeGrid,...
            A_SS,w_SS,rr,NN);
        err = max(abs(mCoefficientsNew(:) - mCoefficients(:)));
        iteration = iteration + 1;
        mCoefficients = dampening * mCoefficients + (1 - dampening) * mCoefficientsNew;

    end

    
    %----------------------------------------------------------------
    % Compute histogram approximation of stationary distribution
    %----------------------------------------------------------------

    %%%
    % Compute policies over histogram grid	
    %%%

    % Compute decision rules along fine grid
    mConditionalExpectation = exp(mCoefficients * vAssetsPolyFine');

    % Compute savings policy
    mAssetsPrimeStar = ((1-ttau)*w_SS*mzGridFine).^(1+1/nnu).*(mConditionalExpectation/ppsi).^(1/nnu) ...
                        +(1+rr)*mAssetsGridFine+chidT_SS-1./mConditionalExpectation;
    mAssetsPrimeFine = max(mAssetsPrimeStar,bbBar);

    % Compute labor and consumption
    % ASSUMES nnu=1!
    mConstr = (mAssetsPrimeFine==bbBar); % Constrained states
    mLaborFine = (1-ttau)*w_SS*mzGridFine.*mConditionalExpectation/ppsi; % Labor supply if not constrained
    aux = -bbBar + (1+rr)*mAssetsGridFine(mConstr) + chidT_SS;
    mLaborFine(mConstr) = (-aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*mzGridFine(mConstr)).^2/ppsi)) ...
                           ./ (2*(1-ttau)*w_SS*mzGridFine(mConstr)); % Labor supply if constrained
    mConsumptionFine = (1-ttau)*w_SS*mzGridFine./(ppsi*mLaborFine); % Consumption

    %%%
    % Compute transition matrix associated with policy rules
    %%%

    % Compute weighting matrices
    [vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGridFine,mAssetsPrimeFine(:));

    % Compute transition matrix for assets over full grid
    mTransitionAbove = zeros(nStateFine,nAssetsFine);
    mTransitionBelow = zeros(nStateFine,nAssetsFine);
    for b = 1:nAssetsFine
        mTransitionBelow(vIndicesBelow == b,b) = vWeightBelow(vIndicesBelow == b);
        mTransitionAbove(vIndicesAbove == b,b) = vWeightAbove(vIndicesAbove == b);
    end
    mAssetsTransition = kron(mTransitionBelow + mTransitionAbove,ones(1,nz));

    % Compute transition matrix for idiosyncratic shocks over full grid
    mzTransitionHistogram = repmat(mzTransition,nAssetsFine);

    % Compute full transition matrix
    mTransition = sparse(mAssetsTransition .* mzTransitionHistogram);

    % Compute invariant histogram by iteration
    errHistogram = Inf;	iterationHistogram = 0;
    vHistogram = ones(nStateFine,1) / nStateFine;
    while errHistogram > 1e-12 && iterationHistogram < maxIterations

        vHistogramNew = mTransition' * vHistogram;
        errHistogram = max(abs(vHistogramNew - vHistogram));
        iterationHistogram = iterationHistogram + 1;
        vHistogram = vHistogramNew;

    end

    % Expand histogram matrix
    mHistogram = reshape(full(vHistogramNew),nz,nAssetsFine);
    
    % Aggregate bond demand and labor supply for consumer type
    bond_demand(iShare) = sum(mHistogram,1) * vAssetsGridFine;
    labor_supply(iShare) = (mzGridFine(:).*mLaborFine(:))' * mHistogram(:);
    
    % Store results
    if nargout > 1
        mHistogramOptional(iShare,:,:) = mHistogram;
        mAssetsPrimeOptional(iShare,:,:) = mAssetsPrimeFine;
        mConsumptionOptional(iShare,:,:) = mConsumptionFine;
        mLaborOptional(iShare,:,:) = mLaborFine;
        mConditionalExpectationOptional(iShare,:,:) = mConditionalExpectation_hist;
    end
    
end

%----------------------------------------------------------------
% Return market clearing residuals
%----------------------------------------------------------------

residual = [vvarthetaB*A_SS*NN + bond_demand'*vShareFraction; % Bonds
            NN - labor_supply'*vShareFraction];    % Labor
