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
nz = var_array{11};
nAssets = var_array{12};
nState = var_array{13};
nStateFine = var_array{14};
nAssetsFine = var_array{15};
nAssetsQuadrature = var_array{16};
nMeasure = var_array{17};
assetsMin = var_array{18};
assetsMax = var_array{19};
vzGrid = var_array{20};
mzGrid = var_array{21};
mzGridFine = var_array{22};
mzPrimeGrid = var_array{23};
mzGridQuadrature = var_array{24};
vAssetsGridFine = var_array{25};
vAssetsGridQuadrature = var_array{26};
vAssetsPoly = var_array{27};
vAssetsPolySquared = var_array{28};
vAssetsPolyFine = var_array{29};
vAssetsPolyQuadrature = var_array{30};
vAssetsPolyBC = var_array{31};
mAssetsGrid = var_array{32};
mAssetsGridFine = var_array{33};
mAssetsGridQuadrature = var_array{34};
vQuadratureWeights = var_array{35};
maxIterations = var_array{36};
tolerance = var_array{37};
dampening = var_array{38};
A_SS = var_array{39};
w_SS = var_array{40};



%----------------------------------------------------------------
% Compute individual decisions
%----------------------------------------------------------------

chidT_SS = (1/eepsilon + vvarthetaT)*A_SS*NN;

% Initialize coefficients using rule of thumb savings rule that sets b'=b and z'=z
% ASSUMES nnu=1!
% aux = rr*mAssetsGrid + chidT_SS;
% mGridInit = log(bbeta * (1+rr) * 2 ./ (aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*mzGrid).^2/ppsi)));
mGridInit = log(mConditionalExpectation);

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

mCoefficientsOptional = mCoefficients;


%----------------------------------------------------------------
% Compute policies over quadrature grid for integration
%----------------------------------------------------------------

% Compute conditional expectation
mConditionalExpectation = exp(mCoefficients * vAssetsPolyQuadrature');

% Compute savings policy
mAssetsPrimeStar = ((1-ttau)*w_SS*mzGridQuadrature).^(1+1/nnu).*(mConditionalExpectation/ppsi).^(1/nnu) ...
                    +(1+rr)*mAssetsGridQuadrature+chidT_SS-1./mConditionalExpectation;
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

% options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter-detailed',...
%     'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);


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
%          objectiveFunction = @(vParametersTilde) parametersResidual(vParametersTilde,reshape(bGridMoments(iz,:,:),nAssetsQuadrature,nMeasure),vQuadratureWeights,nMeasure);
%         [vParameters,normalization] = fminunc(objectiveFunction,zeros(nMeasure,1),options);        

		vParameters = mParameters(iz,2:end)';
		mGridMoments = reshape(bGridMoments(iz,:,:),nAssetsQuadrature,nMeasure);
        dens_nonnormaliz = exp(mGridMoments * vParameters);
		DD = vQuadratureWeights' * (dens_nonnormaliz.*mGridMoments);
		HH = zeros(nMeasure);
        for iH = 1:nAssetsQuadrature
			HH = HH+vQuadratureWeights(iH)*dens_nonnormaliz(iH,:)*(mGridMoments(iH,:)'*mGridMoments(iH,:));
        end
		vParameters = vParameters-HH\DD';
		normalization = vQuadratureWeights' * (exp(mGridMoments * vParameters));
		mParametersNew(iz,:) = [1 / normalization; vParameters];
    end
	
	% Compute new moments and centered moments grid
	mMomentsNew = zeros(nz,nMeasure);
	bGridMomentsNew = zeros(nz,nAssetsQuadrature,nMeasure);
	
	for iz = 1 : nz
		
		% Compute first moment (uncentered)
		mMomentsNew(iz,1) = 0;
		for izTilde = 1 : nz
		
			mMomentsNew(iz,1) = mMomentsNew(iz,1) + (1 - mHat(izTilde,1)) * vzInvariant(izTilde) * mzTransition(...
				izTilde,iz) * mParametersNew(izTilde,1) * vQuadratureWeights' * (mAssetsPrimeQuadrature(izTilde,:)' .* ...
				exp(reshape(bGridMoments(izTilde,:,:),nAssetsQuadrature,nMeasure) * mParametersNew(izTilde,2:nMeasure+1)')) + mHat(izTilde,1) * ...
				vzInvariant(izTilde) * mzTransition(izTilde,iz) * mAssetsPrimeBC(izTilde,1);	
				
		end
		
		mMomentsNew(iz,1) = mMomentsNew(iz,1) / vzInvariant(iz);
		bGridMomentsNew(iz,:,1) = vAssetsGridQuadrature - mMomentsNew(iz,1);
		
		% Compute higher order moments (centered)
		for iMoment = 2 : nMeasure
		
			mMomentsNew(iz,iMoment) = 0;
			
			for izTilde = 1 : nz
			
				mMomentsNew(iz,iMoment) = mMomentsNew(iz,iMoment) + (1 - mHat(izTilde,1)) * vzInvariant(izTilde) * mzTransition(...
					izTilde,iz) * mParametersNew(izTilde,1) * vQuadratureWeights' * (((mAssetsPrimeQuadrature(izTilde,:)' - ...
					mMomentsNew(iz,1)) .^ iMoment) .* exp(reshape(bGridMoments(izTilde,:,:),nAssetsQuadrature,nMeasure) * ...
					mParametersNew(izTilde,2:nMeasure+1)')) + mHat(izTilde,1) * vzInvariant(izTilde) * ...
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
		
			mHatNew(iz,1) = mHatNew(iz,1) + (1 - mHat(izTilde,1)) * vzInvariant(izTilde) * ...
				mzTransition(izTilde,iz) * mParametersNew(izTilde,1) * vQuadratureWeights' * ...
				((mAssetsPrimeQuadrature(izTilde,:)' <= bbBar + 1e-8) .* exp(reshape(bGridMoments(izTilde,:,:),nAssetsQuadrature,nMeasure) * mParametersNew(izTilde,2:nMeasure+1)')) + ...
				mHat(izTilde,1) * vzInvariant(izTilde) * mzTransition(izTilde,iz) * ...
				(mAssetsPrimeBC(izTilde,1) <= bbBar + 1e-8);
				
		end
		
		mHatNew(iz,1) = mHatNew(iz,1) / vzInvariant(iz);
	
    end
	
	%%%
	% Update iteration
	%%%
	
	err = max([max(abs(mMomentsNew(:) - mMoments(:))),max(abs(mHatNew(:) - mHat(:))),max(abs(mParametersNew(:) - mParameters(:)))]);
	iteration = iteration + 1;
	mMoments = mMomentsNew;
	bGridMoments = bGridMomentsNew;
	mHat = mHatNew;
    mParameters = mParametersNew;
	
end

disp(iteration);

%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

% Household ss borrowing
bbNew = (vzInvariant .* (1 - mHat))' * mMoments(:,1) + bbBar * (vzInvariant' * mHat);

% Household ss effective labor supply (=z*n)
mLabor_times_z = (1-ttau)*w_SS*(vzGrid.^2).*mConditionalExpectation/ppsi; % Away from constraint
aux = rr*bbBar + chidT_SS;
mLabor_times_z_BC = (-aux + sqrt(aux.^2 + 4*((1-ttau)*w_SS*vzGrid).^2/ppsi)) ...
                       ./ (2*(1-ttau)*w_SS); % At constraint
NNNew = (vzInvariant.*(1-mHat))'*mLabor_times_z*vQuadratureWeights + (vzInvariant.*mHat)'*mLabor_times_z_BC; % Aggregate labor supply

residual = [vvarthetaB*A_SS*NN + bbNew;
            NN - NNNew];

% Also return optional outputs if requested
if nargout > 2

    mParametersOptional = mParameters;
    mMomentsOptional = mMoments;
	mHatOptional = mHat;
	
end
