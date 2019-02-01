function [residual,mCoefficientsOptional,mParametersOptional,mMomentsOptional,mHatOptional] = ...
	computeMCResidualPolynomials(capital,mMoments,aGridMoments,mHat,var_array)

% Computes residual of market-clearing condition, parametric family to approximate distribution
% 
% Inputs
%   (1) capital: candidate aggregate capital stock
%	(2) mMoments: intial guess of moments of distribution (nEpsilon x nMeasure)
%	(3) aGridMoments: grid of centralized moments for computing PDF, corresponding to mMoments
%		(nEpsilon x nAssetsQuadrature x nMoments)
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
	
% Compute prices
r = aalpha * (capital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (capital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));
var_array{41} = r;
var_array{42} = w;

%----------------------------------------------------------------
% Set error tolerance & max iteration depending on use
%----------------------------------------------------------------

if nargout == 1 
    err1 = tolerance;
    err2 = 1e-4;
    tol2 = 2000;
elseif nargout > 1 
    err1 = 1e-8;
    err2 = 1e-6;
    tol2 = 5000;
else
    error('Check # of inputs & outputs for computeMCResidualPolynomials');
end

%----------------------------------------------------------------
% Compute individual decisions
%----------------------------------------------------------------

if splineOpt == 0	% approximate conditional expectation function using polynomials

	% Initialize coefficients using rule of thumb savings rule
	mGridInit = log(bbeta * (1 + r) * ((w * (mmu * (1 - mEpsilonGrid) + (1 - ttau) * mEpsilonGrid) + ...
		r * mAssetsGrid) .^ (-ssigma)));
	mCoefficients = zeros(nEpsilon,nAssets);
	for iEpsilon = 1:nEpsilon	% interpolate
		vCoefficients = sum(vAssetsPoly' .* (ones(nAssets,1) * mGridInit(iEpsilon,:)),2);
		mCoefficients(iEpsilon,:) = (vCoefficients ./ vAssetsPolySquared)';
	end

	% Iterate
	err = 100; iteration = 1;
	while err > tolerance && iteration <= maxIterations

		mCoefficientsNew = updateCoefficients_polynomials(mCoefficients,var_array);
		err = max(abs(mCoefficientsNew(:) - mCoefficients(:)));
		iteration = iteration + 1;
		mCoefficients = dampening * mCoefficients + (1 - dampening) * mCoefficientsNew;

	end
	
	mCoefficientsOptional = mCoefficients;
	
else	% approximate savings decision using linear splines
	
	% Initialize coefficients
	mAssetsPrime = mAssetsGrid;
	
	% Iterate
	err = 100; iteration = 1;
	while err > tolerance && iteration <= maxIterations
		mAssetsPrimeNew = updateCoefficients_splines(mAssetsPrime);
		err = max(abs(mAssetsPrimeNew(:) - mAssetsPrime(:)));
		iteration = iteration + 1;
		mAssetsPrime = dampening * mAssetsPrime + (1 - dampening) * mAssetsPrimeNew;
	end
	
	mCoefficientsOptional = mAssetsPrime;

end

%----------------------------------------------------------------
% Compute policies over quadrature grid for integration
%----------------------------------------------------------------

if splineOpt == 0

	% Compute conditional expectation
	mConditionalExpectation = exp(mCoefficients * vAssetsPolyQuadrature');

	% Compute savings policy
	mAssetsPrimeStar = w * (mmu * (1 - mEpsilonGridQuadrature) + (1 - ttau) * mEpsilonGridQuadrature) + ...
		(1 + r) * mAssetsGridQuadrature - (mConditionalExpectation .^ (-1 / ssigma));
	mAssetsPrimeQuadrature = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,nAssetsQuadrature));
		
else

	% Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,vAssetsGridQuadrature);
		
	% Linear interpolation
	mAssetsPrimeQuadrature = mAssetsPrime(:,vIndicesBelow) .* repmat(vWeightBelow',nEpsilon,1) + ...
		mAssetsPrime(:,vIndicesAbove) .* repmat(vWeightAbove',nEpsilon,1);
	
end

%----------------------------------------------------------------
% Compute policies at borrowing constraint for integration
%----------------------------------------------------------------

if splineOpt == 0

	% Compute conditional expectation
	mConditionalExpectation = exp(mCoefficients * vAssetsPolyBC');

	% Compute savings policy
	mAssetsPrimeStar = w * (mmu * (1 - vEpsilonGrid) + (1 - ttau) * vEpsilonGrid) + ...
		(1 + r) * assetsMin - (mConditionalExpectation .^ (-1 / ssigma));
	mAssetsPrimeBC = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,1));
		
else

	% Compute weights
	[vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = computeLinearWeights(vAssetsGrid,assetsMin);
		
	% Linear interpolation
	mAssetsPrimeBC = mAssetsPrime(:,vIndicesBelow) .* repmat(vWeightBelow',nEpsilon,1) + ...
		mAssetsPrime(:,vIndicesAbove) .* repmat(vWeightAbove',nEpsilon,1);
	
end

%----------------------------------------------------------------
% Compute stationary distribution from these decision rules
%----------------------------------------------------------------

% Initialize iteration
err = 100; iteration = 1; 
% options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
	% 'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
%{ For older versions of MATLAB:
% options = optimset('LargeScale','off','Display','notify-detailed',...
% 	'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
%}

mParameters = zeros(nEpsilon,nMeasure+1);
% mParameters0 = mParameters;
% Iteration
while err > err2 && iteration <= tol2
	
	%%%
	% Update density away from borrowing constraint
	%%%

    mParametersNew = zeros(nEpsilon,nMeasure+1);

	for iEpsilon = 1 : nEpsilon
		% objectiveFunction = @(vParametersTilde) parametersResidual(vParametersTilde,reshape(aGridMoments(iEpsilon,:,:),nAssetsQuadrature,nMeasure),vQuadratureWeights,nMeasure);
		% [vParameters,normalization] = fminunc(objectiveFunction,zeros(nMeasure,1),options);		
		vParameters = mParameters(iEpsilon,2:end)';
		mGridMoments = reshape(aGridMoments(iEpsilon,:,:),nAssetsQuadrature,nMeasure);
        dens_nonnormaliz = exp(mGridMoments * vParameters);
		DD = vQuadratureWeights' * (dens_nonnormaliz.*mGridMoments);
		HH = zeros(nMeasure);
		for iH = 1:nAssetsQuadrature
			HH = HH+vQuadratureWeights(iH)*dens_nonnormaliz(iH,:)*mGridMoments(iH,:)'*mGridMoments(iH,:);
		end
		vParameters = vParameters-HH\DD';
		normalization = vQuadratureWeights' * dens_nonnormaliz;
		mParametersNew(iEpsilon,:) = [1 / normalization; vParameters];
	end
% 	if mod(iteration,100) == 99
% 		disp(['Parameter convergence (iter' num2str(iteration) '-iter' num2str(iteration-99) '): ' num2str(sum((mParameters(:)-mParameters0(:)).^2))])
% 	end
	
	% Compute new moments and centered moments grid
	mMomentsNew = zeros(nEpsilon,nMeasure);
	aGridMomentsNew = zeros(nEpsilon,nAssetsQuadrature,nMeasure);
	
	for iEpsilon = 1 : nEpsilon
		
		% Compute first moment (uncentered)
		mMomentsNew(iEpsilon,1) = 0;
		for iEpsilonTilde = 1 : nEpsilon
		
			mMomentsNew(iEpsilon,1) = mMomentsNew(iEpsilon,1) + (1 - mHat(iEpsilonTilde,1)) * vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(...
				iEpsilonTilde,iEpsilon) * mParametersNew(iEpsilonTilde,1) * vQuadratureWeights' * (mAssetsPrimeQuadrature(iEpsilonTilde,:)' .* ...
				exp(reshape(aGridMoments(iEpsilonTilde,:,:),nAssetsQuadrature,nMeasure) * mParametersNew(iEpsilonTilde,2:nMeasure+1)')) + mHat(iEpsilonTilde,1) * ...
				vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(iEpsilonTilde,iEpsilon) * mAssetsPrimeBC(iEpsilonTilde,1);	
				
		end
		
		mMomentsNew(iEpsilon,1) = mMomentsNew(iEpsilon,1) / vEpsilonInvariant(iEpsilon);
		aGridMomentsNew(iEpsilon,:,1) = vAssetsGridQuadrature - mMomentsNew(iEpsilon,1);
		
		% Compute higher order moments (centered)
		for iMoment = 2 : nMeasure
		
			mMomentsNew(iEpsilon,iMoment) = 0;
			
			for iEpsilonTilde = 1 : nEpsilon
			
				mMomentsNew(iEpsilon,iMoment) = mMomentsNew(iEpsilon,iMoment) + (1 - mHat(iEpsilonTilde,1)) * vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(...
					iEpsilonTilde,iEpsilon) * mParametersNew(iEpsilonTilde,1) * vQuadratureWeights' * (((mAssetsPrimeQuadrature(iEpsilonTilde,:)' - ...
					mMomentsNew(iEpsilon,1)) .^ iMoment) .* exp(reshape(aGridMoments(iEpsilonTilde,:,:),nAssetsQuadrature,nMeasure) * ...
					mParametersNew(iEpsilonTilde,2:nMeasure+1)')) + mHat(iEpsilonTilde,1) * vEpsilonInvariant(iEpsilonTilde) * ...
					mEpsilonTransition(iEpsilonTilde,iEpsilon) * ((mAssetsPrimeBC(iEpsilonTilde,1) - mMomentsNew(iEpsilon,1)) .^ iMoment);
					
			end
			
			mMomentsNew(iEpsilon,iMoment) = mMomentsNew(iEpsilon,iMoment) / vEpsilonInvariant(iEpsilon);
			aGridMomentsNew(iEpsilon,:,iMoment) = (vAssetsGridQuadrature' - mMomentsNew(iEpsilon,1)) .^ iMoment - ...
				mMomentsNew(iEpsilon,iMoment);
				
		end
		
	end

	%%%
	% Update mass at borrowing constraint
	%%%
	
	mHatNew = zeros(nEpsilon,1);
	
	for iEpsilon = 1 : nEpsilon
	
		for iEpsilonTilde = 1 : nEpsilon
		
			mHatNew(iEpsilon,1) = mHatNew(iEpsilon,1) + (1 - mHat(iEpsilonTilde,1)) * vEpsilonInvariant(iEpsilonTilde) * ...
				mEpsilonTransition(iEpsilonTilde,iEpsilon) * mParametersNew(iEpsilonTilde,1) * vQuadratureWeights' * ...
				((mAssetsPrimeQuadrature(iEpsilonTilde,:)' <= aaBar + 1e-8) .* exp(reshape(aGridMoments(iEpsilonTilde,:,:),nAssetsQuadrature,nMeasure) * mParametersNew(iEpsilonTilde,2:nMeasure+1)')) + ...
				mHat(iEpsilonTilde,1) * vEpsilonInvariant(iEpsilonTilde) * mEpsilonTransition(iEpsilonTilde,iEpsilon) * ...
				(mAssetsPrimeBC(iEpsilonTilde,1) <= aaBar + 1e-8);
				
		end
			
		mHatNew(iEpsilon,1) = mHatNew(iEpsilon,1) / vEpsilonInvariant(iEpsilon);
	
   end
	
	%%%
	% Update iteration
	%%%
	
	err = max([max(abs(mMomentsNew(:) - mMoments(:))),max(abs(mHatNew(:) - mHat(:))),max(abs(mParametersNew(:) - mParameters(:)))]);
	iteration = iteration + 1;
	mMoments = mMomentsNew;
	aGridMoments = aGridMomentsNew;
	mHat = mHatNew;
    mParameters = mParametersNew;
	
end

disp(iteration);

%----------------------------------------------------------------
% Return market clearing residual
%----------------------------------------------------------------

capitalNew = (vEpsilonInvariant .* (1 - mHat))' * mMoments(:,1) + aaBar * (vEpsilonInvariant .* mHat)' * ones(nEpsilon,1);
residual = capital - capitalNew;

% Also return optional outputs if requested
if nargout > 2

    mParametersOptional = mParameters;
    mMomentsOptional = mMoments;
	mHatOptional = mHat;
	
end
