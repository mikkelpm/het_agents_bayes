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
var_array{11} = nz;
var_array{12} = nAssets;
var_array{13} = nState;
var_array{14} = nStateFine;
var_array{15} = nAssetsFine;
var_array{16} = nAssetsQuadrature;
var_array{17} = nMeasure;
var_array{18} = assetsMin;
var_array{19} = assetsMax;
var_array{20} = vzGrid;
var_array{21} = mzGrid;
var_array{22} = mzGridFine;
var_array{23} = mzPrimeGrid;
var_array{24} = mzGridQuadrature;
var_array{25} = vAssetsGridFine;
var_array{26} = vAssetsGridQuadrature;
var_array{27} = vAssetsPoly;
var_array{28} = vAssetsPolySquared;
var_array{29} = vAssetsPolyFine;
var_array{30} = vAssetsPolyQuadrature;
var_array{31} = vAssetsPolyBC;
var_array{32} = mAssetsGrid;
var_array{33} = mAssetsGridFine;
var_array{34} = mAssetsGridQuadrature;
var_array{35} = vQuadratureWeights;
var_array{36} = maxIterations;
var_array{37} = tolerance;
var_array{38} = dampening;
var_array{39} = A_SS;
var_array{40} = w_SS;


%----------------------------------------------------------------
% Compute initial guess of market-clearing r and N using
% histogram approximation of distribution, from Young (2010)
%----------------------------------------------------------------

t0 = tic;
fprintf('Computing initial guess from histogram...\n')

% Solve for market clearing capital stock
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
[~, mHistogram, mAssets, mConsumption, mLabor, mConditionalExpectation] = computeMCResidualHistogram(r_SS_hist,N_SS_hist,var_array);

% Compute moments from histogram
mMomentsHistogram = zeros(nz,nMeasure);
bGridMoments = zeros(nz,nAssetsQuadrature,nMeasure); % grid for computing PDF

for iz = 1 : nz
	
	% First moment (uncentered)
	mMomentsHistogram(iz,1) = sum(vAssetsGridFine' .* (mHistogram(iz,:) ./ ...
		sum(mHistogram(iz,:))));
	bGridMoments(iz,:,1) = vAssetsGridQuadrature - mMomentsHistogram(iz,1);
		
	% Higher order moments (centered)
	for iMoment = 2 : nMeasure
		mMomentsHistogram(iz,iMoment) = sum(((vAssetsGridFine' - mMomentsHistogram(iz,1)) .^ iMoment) .* ...
			(mHistogram(iz,:) ./ sum(mHistogram(iz,:))));
		bGridMoments(iz,:,iMoment) = (vAssetsGridQuadrature' - mMomentsHistogram(iz,1)) .^ ...
			iMoment - mMomentsHistogram(iz,iMoment);
	end	
	
end

% Mass at borrowing constraint
mHatHistogram = mHistogram(:,1)./sum(mHistogram,2);

%----------------------------------------------------------------
% Compute market-clearing capital stock from parametric family
%----------------------------------------------------------------

t0 = tic; 
fprintf('Compute steady state from parametric family...\n')

% Solve for market clearing capital stock
f = @(x) computeMCResidualPolynomials(x(1),x(2),mMomentsHistogram,bGridMoments,mHatHistogram,mConditionalExpectation,var_array);
options = optimoptions('fsolve','Display',displayOpt,'TolFun',tolerance_SS_root);
if sum(abs(f([r_SS_hist,N_SS_hist]))) > tolerance_SS_root
    [x_solve,err,exitflag] = fsolve(f,[r_SS_hist,N_SS_hist],options);
    r_SS = x_solve(1);
    N_SS = x_solve(2);
    % Return error if market clearing not solved
    if exitflag < 1
        check = 1;
        return; 
    end	
else    
	r_SS = r_SS_hist;
    N_SS = N_SS_hist;
end

if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end


%----------------------------------------------------------------
% Compute other objects from steady state
%----------------------------------------------------------------

[~,mCoefficients,mParameters,mMoments,mHat] = ...
    computeMCResidualPolynomials(r_SS,N_SS,mMomentsHistogram,bGridMoments,mHatHistogram,mConditionalExpectation,var_array);

Y_SS = A_SS*N_SS;
d_SS = Y_SS/eepsilon;
B_SS = vvarthetaB*Y_SS;
T_SS = vvarthetaT*Y_SS;
G_SS = r_SS*B_SS + ttau*w_SS*N_SS - T_SS;

