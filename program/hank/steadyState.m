% Computes and analyzes steady state with no aggregate shocks
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

profile on

oldFolder = cd('./Auxiliary Functions');

%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Compute Steady State
%----------------------------------------------------------------

% Solve for steady state capital stock, distribution, and decision rules
coreSteadyState;
profile viewer;
return;
% profsave(profile('info'),['profile_results_' datestr(now,'yyyymmdd')]);

% Compute decision rules along fine grid for analysis
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
[~,mHistogram,mAssetsPrime,mConsumption] = computeMCResidualHistogram(aggregateCapital,var_array);

% Compute density along fine grid
mDistributionFine = zeros(nEpsilon,nAssetsFine);
for iEpsilon = 1 : nEpsilon

	% First moment (uncentered)
	mGridMoments = zeros(nAssetsFine,nMeasure);
	mGridMoments(:,1) = (vAssetsGridFine - mMoments(iEpsilon,1));
	
	% Higher order moments (centered)
	for iMoment = 2 : nMeasure
		mGridMoments(:,iMoment) = (vAssetsGridFine - mMoments(iEpsilon,1)) .^ iMoment - ...
			mMoments(iEpsilon,iMoment);
	end
	
	% Compute density away from borrowing constraint
	mDistributionFine(iEpsilon,:) = mParameters(iEpsilon,1) * exp(mGridMoments * ...
		mParameters(iEpsilon,2:nMeasure+1)');
		
	% Mass at borrowing constraint
	%mDistributionFine(iEpsilon,1) = mHat(iEpsilon,1); % Commented out for now; need fine quadrature grid to capture correctly
		
end

%----------------------------------------------------------------
% Plot results 
%----------------------------------------------------------------

% Savings function
figure
hold on
plot(vAssetsGridFine,mAssetsPrime(1,:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vAssetsGridFine,mAssetsPrime(2,:),'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vAssetsGridFine,vAssetsGridFine,'k--','linewidth',1)
xlabel('Assets, $a$','interpreter','latex')
ylabel('Savings, $s(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Savings Decision Rule')
legend('Unemployed','Employed','location','southeast')
grid on
set(gcf,'color','w')
hold off

% Consumption function
figure
hold on
plot(vAssetsGridFine,mConsumption(1,:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vAssetsGridFine,mConsumption(2,:),'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Assets, $a$','interpreter','latex')
ylabel('Consumption, $c(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Consumption Decision rule')
legend('Unemployed','Employed','location','southeast')
grid on
set(gcf,'color','w')
hold off

% Distribution of unemployed households
figure
hold on
plot(vAssetsGridFine,mHistogram(1,:) / sum(mHistogram(1,:)),...
	'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vAssetsGridFine,mDistributionFine(1,:) ./ sum(mDistributionFine(1,:)),...
	'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Invariant Distribution of Households (Unemployed)')
legend('Histogram','Parametric Family','location','northeast')
grid on
set(gcf,'color','w')
hold off

% Distribution of employed households
figure
hold on
plot(vAssetsGridFine,mHistogram(2,:) / sum(mHistogram(2,:)),...
	'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vAssetsGridFine,mDistributionFine(2,:) ./ sum(mDistributionFine(2,:)),...
	'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Invariant Distribution of Households (Employed)')
legend('Histogram','Parametric Family','location','northeast')
grid on
set(gcf,'color','w')
hold off

cd(oldFolder)