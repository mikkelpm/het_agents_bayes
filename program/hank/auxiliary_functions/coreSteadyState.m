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
computePolynomials;

% variables in array to pass on
var_array{1} = bbeta;
var_array{2} = ppsi;
var_array{3} = nnu;
var_array{4} = bbBar;
var_array{5} = eepsilon;
var_array{6} = ttau;
var_array{7} = vvarthetaB;
var_array{8} = vvarthetaT;
var_array{9} = mzTransition;
var_array{10} = vzInvariant;
var_array{11} = vShareGrid;
var_array{12} = vShareFraction;
var_array{13} = nz;
var_array{14} = nAssets;
var_array{15} = nState;
var_array{16} = nStateFine;
var_array{17} = nAssetsFine;
var_array{18} = nAssetsQuadrature;
var_array{19} = nMeasure;
var_array{20} = nShare;
var_array{21} = assetsMin;
var_array{22} = assetsMax;
var_array{23} = vzGrid;
var_array{24} = mzGrid;
var_array{25} = mzGridFine;
var_array{26} = mzPrimeGrid;
var_array{27} = mzGridQuadrature;
var_array{28} = vAssetsGridFine;
var_array{29} = vAssetsGridQuadrature;
var_array{30} = vAssetsPoly;
var_array{31} = vAssetsPolySquared;
var_array{32} = vAssetsPolyFine;
var_array{33} = vAssetsPolyQuadrature;
var_array{34} = vAssetsPolyBC;
var_array{35} = mAssetsGrid;
var_array{36} = mAssetsGridFine;
var_array{37} = mAssetsGridQuadrature;
var_array{38} = vQuadratureWeights;
var_array{39} = maxIterations;
var_array{40} = tolerance;
var_array{41} = dampening;
var_array{42} = numNewton;
var_array{43} = A_SS;
var_array{44} = w_SS;


%----------------------------------------------------------------
% Compute initial guess of market-clearing r and N using
% histogram approximation of distribution, from Young (2010)
%----------------------------------------------------------------

t0 = tic;
fprintf('Computing initial guess from histogram...\n')

% Solve for market clearing r and N
f = @(x) computeMCResidualHistogram(x(1),x(2),var_array);
options = optimoptions('fsolve','Display',displayOpt,...
                       'TolFun',tolerance_SS_root,'TypicalX',[0.01; 0.1]);
[x_solve_hist,err,exitflag] = fsolve(f,[r_RepSS N_RepSS],options);

% Return exitflag if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end	

r_SS_hist = x_solve_hist(1);
N_SS_hist = x_solve_hist(2);
if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%----------------------------------------------------------------
% Compute moments of histogram to use as initial guess for parametric family
%----------------------------------------------------------------

% Compute histogram
[~, mHistogram, mAssetsPrime_hist, mConsumption_hist, mLabor_hist, mConditionalExpectation_hist] ...
    = computeMCResidualHistogram(r_SS_hist,N_SS_hist,var_array);

% Compute moments from histogram
mMomentsHistogram = zeros(nz,nMeasure,nShare);
bGridMoments = zeros(nz,nAssetsQuadrature,nMeasure,nShare); % grid for computing PDF

for iShare = 1:nShare
    for iz = 1 : nz
        
        % First moment (uncentered)
        mMomentsHistogram(iz,1,iShare) = mHistogram(iz,:,iShare) * vAssetsGridFine / ...
            sum(mHistogram(iz,:,iShare));
        bGridMoments(iz,:,1,iShare) = vAssetsGridQuadrature - mMomentsHistogram(iz,1,iShare);
        
        % Higher order moments (centered)
        for iMoment = 2 : nMeasure
            mMomentsHistogram(iShare,iz,iMoment) = mHistogram(iz,:,iShare) *...
			((vAssetsGridFine - mMomentsHistogram(iz,1,iShare)) .^ iMoment) / sum(mHistogram(iz,:,iShare));
            bGridMoments(iz,:,iMoment,iShare) = (vAssetsGridQuadrature' - mMomentsHistogram(iz,1,iShare)) .^ ...
                iMoment - mMomentsHistogram(iz,iMoment,iShare);
        end
        
    end
end

% Mass at borrowing constraint
mHat_hist = mHistogram(:,:,1)./sum(mHistogram,3);

%----------------------------------------------------------------
% Compute market-clearing capital stock from parametric family
%----------------------------------------------------------------

t0 = tic; 
fprintf('Compute steady state from parametric family...\n')

% Solve for market clearing r and N
f = @(x) computeMCResidualPolynomials(x(1),x(2),mMomentsHistogram,bGridMoments,mHat_hist,mConditionalExpectation_hist,var_array);
[x_solve,err,exitflag] = fsolve(f,[r_SS_hist,N_SS_hist],options);
r_SS = x_solve(1);
N_SS = x_solve(2);
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
    computeMCResidualPolynomials(r_SS,N_SS,mMomentsHistogram,bGridMoments,mHat_hist,mConditionalExpectation_hist,var_array);

Y_SS = A_SS*N_SS;
d_SS = Y_SS/eepsilon;
B_SS = vvarthetaB*Y_SS;
T_SS = vvarthetaT*Y_SS;
G_SS = r_SS*B_SS + ttau*w_SS*N_SS - T_SS;

