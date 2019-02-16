function [residual,mCoefficientsOptional,mParametersOptional,mMomentsOptional,mHatOptional] = ...
	computeMCResidualPolynomials(rr,NN,mMoments,bGridMoments,mHat,mConditionalExpectation,var_array)

% Computes residual of market-clearing condition, parametric family to approximate distribution
% 
% Inputs
%   (1) capital: candidate aggregate capital stock
%	(2) mMoments: intial guess of moments of distribution (nz x nMeasure)
%	(3) aGridMoments: grid of centralized moments for computing PDF, corresponding to mMoments
%		(nz x nAssetsQuadrature x nMoments)
%	(4) mHat: initial guess of mass at borrowing constraint
%
% Outputs
%   (1) residual: residual of market clearing condition
%   (2) (optional) mCoefficientsOptional: coefficients on individual decisions (splines or polynomials)
%   (3) (optional) mParametersOptional: parameters of density away from borrowing constraint
%   (4) (optional) mMomentsOptional: moments of density away from borrowing constraint
%	(5) (optional) mHatOptional: mass at borrowing constraint
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


mCoefficientsOptional = nan(nz,nAssets,nShare);
mParametersOptional = nan(nz,nMeasure+1,nShare);
mMomentsOptional = nan(nz,nMeasure,nShare);
mHatOptional = nan(nz,nShare);
bond_demand = nan(nShare,1);
labor_supply = nan(nShare,1);

for iShare = 1:nShare

    the_mMoments = mMoments(:,:,iShare);
    the_bGridMoments = bGridMoments(:,:,:,iShare);
    the_mHat = mHat(:,iShare);
    
    %----------------------------------------------------------------
    % Compute individual decisions
    %----------------------------------------------------------------

    chidT_SS = (1/eepsilon*vShareGrid(iShare) + vvarthetaT)*A_SS*NN; % chi*d+T

    % Initialize coefficients using rule of thumb savings rule that sets b'=b and z'=z
    % ASSUMES nnu=1!
    % aux = rr*mAssetsGrid + chidT_SS;
    % mGridInit = log(bbeta * (1+rr) * 2 ./ (aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*mzGrid).^2/ppsi)));

    % Initialize using histogram see
    mGridInit = log(mConditionalExpectation(:,:,iShare));

    mCoefficients = zeros(nz,nAssets);
    for iz = 1:nz	% interpolate
        vCoefficients = vAssetsPoly' * mGridInit(iz,:)';
        mCoefficients(iz,:) = (vCoefficients ./ vAssetsPolySquared)';
    end

    % Iterate
    err = Inf; iteration = 1;
    while err > tolerance && iteration <= maxIterations

        mCoefficientsNew = updateCoefficients_polynomials(mCoefficients,...
            bbeta,ppsi,nnu,bbBar,eepsilon,ttau,vvarthetaT,mzTransition,nz,nAssets,nState,assetsMin,assetsMax,...
            mzGrid,mAssetsGrid,vAssetsPoly,vAssetsPolySquared,mzPrimeGrid,...
            A_SS,w_SS,rr,NN);
        err = max(abs(mCoefficientsNew(:) - mCoefficients(:)));
        iteration = iteration + 1;
        mCoefficients = dampening * mCoefficients + (1 - dampening) * mCoefficientsNew;

    end


    %----------------------------------------------------------------
    % Compute policies over quadrature grid for integration
    %----------------------------------------------------------------

    % Compute conditional expectation
    the_mConditionalExpectation = exp(mCoefficients * vAssetsPolyQuadrature');

    % Compute savings policy
    mAssetsPrimeStar = ((1-ttau)*w_SS*mzGridQuadrature).^(1+1/nnu).*(the_mConditionalExpectation/ppsi).^(1/nnu) ...
                        +(1+rr)*mAssetsGridQuadrature+chidT_SS-1./the_mConditionalExpectation;
    mAssetsPrimeQuadrature = max(mAssetsPrimeStar,bbBar);


    %----------------------------------------------------------------
    % Compute policies at borrowing constraint for integration
    %----------------------------------------------------------------

    % Compute conditional expectation
    mConditionalExpectationBC = exp(mCoefficients * vAssetsPolyBC');

    % Compute savings policy
    mAssetsPrimeStarBC = ((1-ttau)*w_SS*vzGrid).^(1+1/nnu).*(mConditionalExpectationBC/ppsi).^(1/nnu) ...
                        +(1+rr)*assetsMin+chidT_SS-1./mConditionalExpectationBC;
    mAssetsPrimeBC = max(mAssetsPrimeStarBC,bbBar);


    %----------------------------------------------------------------
    % Compute stationary distribution from these decision rules
    %----------------------------------------------------------------

    % Initialize iteration
    err = Inf; iteration = 1; 
    mParameters = zeros(nz,nMeasure+1);

    % Iteration
    while err > tolerance && iteration <= maxIterations

        %%%
        % Update density away from borrowing constraint
        %%%

        mParametersNew = zeros(nz,nMeasure+1);

        for iz = 1 : nz

            vParameters = (1-dampening)*mParameters(iz,2:end)'; %zeros(nMeasure,1);
            mGridMoments = reshape(the_bGridMoments(iz,:,:),nAssetsQuadrature,nMeasure);

            % Take Newton steps to minimize density
            for istep=1:numNewton
                aux = mGridMoments * vParameters;
                dens_nonnormaliz = exp(aux-max(aux(:))); % Non-normalized density
                DD = vQuadratureWeights' * (dens_nonnormaliz.*mGridMoments);
                HH = 0.1*eye(nMeasure); %zeros(nMeasure);
                for iH = 1:nAssetsQuadrature
                    HH = HH+vQuadratureWeights(iH)*dens_nonnormaliz(iH)*(mGridMoments(iH,:)'*mGridMoments(iH,:));
                end
                vParameters = vParameters-HH\DD';
            end

            normalization = vQuadratureWeights' * exp(mGridMoments * vParameters); % New normalization constant
            mParametersNew(iz,:) = [1 / normalization; vParameters]; % New parameters

        end

        % Compute new moments and centered moments grid
        mMomentsNew = zeros(nz,nMeasure);
        bGridMomentsNew = zeros(nz,nAssetsQuadrature,nMeasure);

        for iz = 1 : nz

            % Compute first moment (uncentered)
            mMomentsNew(iz,1) = 0;
            for izTilde = 1 : nz

                mMomentsNew(iz,1) = mMomentsNew(iz,1) + (1 - the_mHat(izTilde,1)) * vzInvariant(izTilde) * mzTransition(...
                    izTilde,iz) * mParametersNew(izTilde,1) * vQuadratureWeights' * (mAssetsPrimeQuadrature(izTilde,:)' .* ...
                    exp(reshape(the_bGridMoments(izTilde,:,:),nAssetsQuadrature,nMeasure) * mParametersNew(izTilde,2:nMeasure+1)')) + the_mHat(izTilde,1) * ...
                    vzInvariant(izTilde) * mzTransition(izTilde,iz) * mAssetsPrimeBC(izTilde,1);	

            end

            mMomentsNew(iz,1) = mMomentsNew(iz,1) / vzInvariant(iz);
            bGridMomentsNew(iz,:,1) = vAssetsGridQuadrature - mMomentsNew(iz,1);

            % Compute higher order moments (centered)
            for iMoment = 2 : nMeasure

                mMomentsNew(iz,iMoment) = 0;

                for izTilde = 1 : nz

                    mMomentsNew(iz,iMoment) = mMomentsNew(iz,iMoment) + (1 - the_mHat(izTilde,1)) * vzInvariant(izTilde) * mzTransition(...
                        izTilde,iz) * mParametersNew(izTilde,1) * vQuadratureWeights' * (((mAssetsPrimeQuadrature(izTilde,:)' - ...
                        mMomentsNew(iz,1)) .^ iMoment) .* exp(reshape(the_bGridMoments(izTilde,:,:),nAssetsQuadrature,nMeasure) * ...
                        mParametersNew(izTilde,2:nMeasure+1)')) + the_mHat(izTilde,1) * vzInvariant(izTilde) * ...
                        mzTransition(izTilde,iz) * ((mAssetsPrimeBC(izTilde,1) - mMomentsNew(iz,1)) .^ iMoment);

                end

                mMomentsNew(iz,iMoment) = mMomentsNew(iz,iMoment) / vzInvariant(iz);
                bGridMomentsNew(iz,:,iMoment) = (vAssetsGridQuadrature' - mMomentsNew(iz,1)) .^ iMoment - ...
                    mMomentsNew(iz,iMoment);

            end

        end

        %%%
        % Update mass at borrowing constraint
        %%%

        mHatNew = zeros(nz,1);

        for iz = 1 : nz

            for izTilde = 1 : nz

                mHatNew(iz,1) = mHatNew(iz,1) + (1 - the_mHat(izTilde,1)) * vzInvariant(izTilde) * ...
                    mzTransition(izTilde,iz) * mParametersNew(izTilde,1) * vQuadratureWeights' * ...
                    ((mAssetsPrimeQuadrature(izTilde,:)' <= bbBar + 1e-8) .* exp(reshape(the_bGridMoments(izTilde,:,:),nAssetsQuadrature,nMeasure) * mParametersNew(izTilde,2:nMeasure+1)')) + ...
                    the_mHat(izTilde,1) * vzInvariant(izTilde) * mzTransition(izTilde,iz) * ...
                    (mAssetsPrimeBC(izTilde,1) <= bbBar + 1e-8);

            end

            mHatNew(iz,1) = mHatNew(iz,1) / vzInvariant(iz);

        end

        %%%
        % Update iteration
        %%%

        err = max([max(abs(mMomentsNew(:) - the_mMoments(:))),max(abs(mHatNew(:) - the_mHat(:))),max(abs(mParametersNew(:) - mParameters(:)))]);
        iteration = iteration + 1;
        the_mMoments = mMomentsNew;
        the_bGridMoments = bGridMomentsNew;
        the_mHat = mHatNew;
        mParameters = mParametersNew;

    end

    if iteration>maxIterations
        warning('Did not converge');
    end
    
    % Compute aggregate bond demand and labor supply
    
    % Household ss borrowing
    bond_demand(iShare) = (vzInvariant .* (1 - the_mHat))' * the_mMoments(:,1) + bbBar * (vzInvariant' * the_mHat);

    % Household ss effective labor supply (=z*n) for b>bbBar
    mConstr = (mAssetsPrimeQuadrature==bbBar);
    mLabor_times_z = (1-ttau)*w_SS*vzGrid.*the_mConditionalExpectation/ppsi;
    aux = -bbBar + (1+rr)*mAssetsPrimeQuadrature(mConstr) + chidT_SS;
    mLabor_times_z(mConstr) = (-aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*mzGridQuadrature(mConstr)).^2/ppsi)) ...
                              ./ (2*(1-ttau)*w_SS); % If savings constrained

    % Household ss effective labor supply for b=bbBar
    mConstrBC = (mAssetsPrimeBC==bbBar);
    mLabor_times_zBC = (1-ttau)*w_SS*vzGrid.*mConditionalExpectationBC/ppsi;
    aux = -bbBar + (1+rr)*mAssetsPrimeBC(mConstrBC) + chidT_SS;
    mLabor_times_zBC(mConstrBC) = (-aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*vzGrid(mConstrBC)).^2/ppsi)) ...
                                  ./ (2*(1-ttau)*w_SS); % If savings constrained

    % Log asset density away from constraint
    logdens = zeros(nz,nAssetsQuadrature);
    for iMoment=1:nMeasure
        logdens = logdens + the_bGridMoments(:,:,iMoment).*mParameters(:,iMoment+1);
    end
    
    % Aggregate effective labor supply
    labor_supply(iShare) = (vzInvariant.*(1-the_mHat).*mParameters(:,1))'*(mLabor_times_z.*exp(logdens))*vQuadratureWeights ...
        + (vzInvariant.*the_mHat)'*mLabor_times_zBC;
    
    % Store results
    if nargout > 1
        mCoefficientsOptional(:,:,iShare) = mCoefficients;
        mParametersOptional(:,:,iShare) = mParameters;
        mMomentsOptional(:,:,iShare) = the_mMoments;
        mHatOptional(:,iShare) = the_mHat;
    end
    
end


%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

residual = [vvarthetaB*A_SS*NN + bond_demand'*vShareFraction; % Bonds
            NN - labor_supply'*vShareFraction];               % Labor

