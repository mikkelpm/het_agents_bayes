clear all;


%% Settings

bbetas = [.95 .96];%linspace(0.9,0.99,10);     % beta values to loop over
T = 200;                            % Number of periods of simulated data
num_burnin_periods = 100;           % Number of burn-in periods for simulations
num_smooth_draws = 25;              % Number of draws from the smoothing distribution (for each beta)
rng_seed = 20180116;                % Random number generator seed for initial simulation


%% Set economic parameters 

global bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment uDuration ...
	mmu rrhoTFP ssigmaTFP;
	
% Preferences
bbeta = .96;										% discount factor (annual calibration)
ssigma = 1;											% coefficient of relative risk aversion
aaBar = 0;											% borrowing constraint

% Technology
aalpha = .36;										% capital share
ddelta = .1;										% depreciation rate (annual calibration)

% Idiosyncratic shocks
vEpsilonGrid = [0;1];
aggEmployment = .93;
uDuration = 1;

% Unemployment benefits
mmu = .15;

% Aggregate Shocks
rrhoTFP = .859;										
ssigmaTFP = .014;


%% Set approximation parameters

global nEpsilon nAssets nAssetsFine nAssetsQuadrature ...
	nMeasure maxIterations tolerance dampening splineOpt displayOpt;

% Whether approximating decision rule with splines or polynomials
splineOpt = 0;	% if splineOpt = 1, use splines to approximate savings policy; if splineOpt = 0, use polynomials
				% to approximate conditional expectation function

% Whether to print out results from steady state computation
displayOpt = 'off';       % 'iter-detailed' or 'off'

% Order of approximation
nEpsilon = 2;
nAssets = 25; % number of gridpoints in spline approximation or polynomials in polynomial approximation

% Finer grid for analyzing policy functions
nAssetsFine = 100;

% Approximation of distribution
nMeasure = 3;
nAssetsQuadrature = 8;

% Iteration on individual decisions
maxIterations = 2e4;
tolerance = 1e-5;
dampening = .95;


%% Save parameters

cd('./Auxiliary Functions');
saveParameters;


%% Initial Dynare run

dynare firstOrderDynamics_polynomials noclearall;                   % Run Dynare once to process model file


%% Simulate data

% Simulate
set_dynare_seed(rng_seed);                                          % Seed RNG
sim_struct = simulate_model(T,num_burnin_periods,M_,oo_,options_);  % Simulate data

% Store simulated data
save('simul.mat', '-struct', 'sim_struct');                         % Save simulated data

%% Smoothing and likelihood

% Determine variables to smooth
smooth_vars_aux1 = cell(nEpsilon,nMeasure);
smooth_vars_aux2 = cell(nEpsilon,nMeasure);
for i_Measure = 1:nMeasure
    for i_Epsilon = 1:nEpsilon
        smooth_vars_aux1{i_Epsilon,i_Measure} = ['moment_' num2str(i_Epsilon) '_' num2str(i_Measure)];
        smooth_vars_aux2{i_Epsilon,i_Measure} = ['measureCoefficient_' num2str(i_Epsilon) '_' num2str(i_Measure)];
    end
end

smooth_vars = char([{'w'; 'r'; 'mHat_1' ; 'mHat_2'}; smooth_vars_aux1(:); smooth_vars_aux2(:)]);%char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones


% Loop over beta values
loglikes = zeros(size(bbetas));
smooth_draws = cell(length(bbetas),num_smooth_draws);

% Heterogeneous agent
ttau = mmu*(1-aggEmployment)/aggEmployment;
loglikes2 = zeros(size(bbetas));
N_hh = 1000;
n_obs_hh = 2; % [labor, income]

for i_beta=[2 1]%1:length(bbetas) % For each alpha...

    fprintf('%s%6.4f\n', 'beta=', bbetas(i_beta));
    bbeta = bbetas(i_beta);                          % Set beta
    saveParameters;
    setDynareParameters;

    % Durbin-Koopman simulation smoother
    [the_ll, ~, the_smooth_draws] = simulation_smoother('simul.mat', smooth_vars, num_smooth_draws, num_burnin_periods, M_, oo_, options_);
    loglikes(i_beta) = the_ll;                       % Store log likelihood
    smooth_draws(i_beta,:) = the_smooth_draws;       % Store smoothing draws

    %% Simulate individual labor and income
    
    if i_beta == 2
        % assign individuals to different smooth draws
        simul_data_hh_grpid = nan(2,N_hh,T);
        simul_data_hh_grpid(1,:,:) = 1; % labor = 1
        N_unemp_hh = round((1-aggEmployment)*N_hh);
        simul_data_hh_grpid(1,1:N_unemp_hh,:) = 0; % labor = 0
        simul_data_hh_grpid(2,:,:) = repmat(1:num_smooth_draws,1,N_hh/num_smooth_draws,T); % id of which smooth draws
        simul_data_hh_grpid(2,floor(N_unemp_hh/num_smooth_draws)*num_smooth_draws+(1:num_smooth_draws),:)...
            =  reshape(randsample(num_smooth_draws,num_smooth_draws*T,true),num_smooth_draws,T);
        
        % draw individual incomes
        simul_data_hh = nan(n_obs_hh,N_hh,T);
        simul_data_hh(1,:,:) = simul_data_hh_grpid(1,:,:);
        for t = 1:T
            for eepsilon = 0:1
                for i_grp = 1:num_smooth_draws
                    ix = find(simul_data_hh_grpid(1,:,t)==eepsilon & simul_data_hh_grpid(2,:,t)==i_grp);
                    
                    % prepare the parameters
                    mHat = eval(['smooth_draws{' num2str(i_beta) ',' num2str(i_grp) '}.mHat_' num2str(eepsilon+1) '(t)']);
                    moment = nan(1,nMeasure);
                    measureCoefficient = nan(1,nMeasure);
                    for i_Measure = 1:nMeasure
                        moment(i_Measure) = eval(['smooth_draws{' num2str(i_beta) ',' num2str(i_grp) '}.moment_' num2str(eepsilon+1) '_' num2str(i_Measure) '(t)']);
                        measureCoefficient(i_Measure) = eval(['smooth_draws{' num2str(i_beta) ',' num2str(i_grp) '}.measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure) '(t)']);
                    end
                    
                    % draw asset
                    n_ix = length(ix);
                    asset_aux = nan(1,n_ix);
                    ix0 = rand(1,n_ix)<mHat;
                    ix1 = setdiff(1:n_ix,ix0);
                    
                    % constrained assets
                    asset_aux(ix0) = aaBar;
                    
                    % unconstrained assets
                    mGridMoments = zeros(nAssetsQuadrature,nMeasure);
                    mGridMoments(:,1) = vAssetsGridQuadrature-moment(1);
                    for iMoment = 2:nMeasure
                        mGridMoments(:,iMoment) = (vAssetsGridQuadrature-moment(1)).^iMoment-moment(iMoment);
                    end
                    normalization = parametersResidual(measureCoefficient',mGridMoments);
                    
                    n_ix1 = length(ix1);
                    uu = rand(1,n_ix1);
                    ccdf = @(a) integral(@(x) exp((x-moment(1))*measureCoefficient(1)...
                        +((x-moment(1)).^2-moment(2))*measureCoefficient(2)...
                        +(x-moment(1)).^3-moment(3)*measureCoefficient(3))/normalization,-inf,a);
                    for i_ix1 = 1:n_ix1
                        asset_aux(ix1(i_ix1)) = fzero(@(a) ccdf(a)-uu(i_ix1),moment(1));
                    end
                    
                    if eepsilon == 0
                        simul_data_hh(2,ix,t) = smooth_draws{i_beta,i_grp}.w(t)*mmu+(1+smooth_draws{i_beta,i_grp}.r(t))*asset_aux;
                    else
                        simul_data_hh(2,ix,t) = smooth_draws{i_beta,i_grp}.w(t)*(1-ttau)+(1+smooth_draws{i_beta,i_grp}.r(t))*asset_aux;
                    end
                end
            end
        end
        save('simul_data_hh_beta.mat','simul_data_hh_grpid','simul_data_hh')   
    end
    %% Likelihood
   
    loglikes_hh = nan(T,N_hh);
    for t = 1:T
        for eepsilon = 0:1
            ix = find(simul_data_hh_grpid(1,:,t)==eepsilon);
            prob_aux = nan(num_smooth_draws,length(ix));
            for i_grp = 1:num_smooth_draws
                
                % assets
                if eepsilon == 0
                    asset_aux = (simul_data_hh(2,ix,t)-smooth_draws{i_beta,i_grp}.w(t)*mmu)/(1+smooth_draws{i_beta,i_grp}.r(t));
                else
                    asset_aux = (simul_data_hh(2,ix,t)-smooth_draws{i_beta,i_grp}.w(t)*(1-ttau))/(1+smooth_draws{i_beta,i_grp}.r(t));
                end
                
                % probability
                % prepare the parameters
                ix0 = abs(asset_aux-aaBar)<1e-6;
                ix1 = abs(asset_aux-aaBar)>=1e-6;
                mHat = eval(['smooth_draws{' num2str(i_beta) ',' num2str(i_grp) '}.mHat_' num2str(eepsilon+1) '(t)']);
                moment = nan(1,nMeasure);
                measureCoefficient = nan(1,nMeasure);
                for i_Measure = 1:nMeasure
                    moment(i_Measure) = eval(['smooth_draws{' num2str(i_beta) ',' num2str(i_grp) '}.moment_' num2str(eepsilon+1) '_' num2str(i_Measure) '(t)']);
                    measureCoefficient(i_Measure) = eval(['smooth_draws{' num2str(i_beta) ',' num2str(i_grp) '}.measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure) '(t)']);
                end
                
                % constrained assets
                prob_aux(i_grp,ix0) = mHat;
                
                % unconstrained assets
                mGridMoments = zeros(nAssetsQuadrature,nMeasure);
                mGridMoments(:,1) = vAssetsGridQuadrature-moment(1);
                for iMoment = 2:nMeasure
                    mGridMoments(:,iMoment) = (vAssetsGridQuadrature-moment(1)).^iMoment-moment(iMoment);
                end
                normalization = parametersResidual(measureCoefficient',mGridMoments);
                
                mGridMoments = zeros(length(ix1),nMeasure);
                mGridMoments(:,1) = asset_aux(ix1)'-moment(1);
                for iMoment = 2:nMeasure
                    mGridMoments(:,iMoment) = (asset_aux(ix1)'-moment(1)).^iMoment-moment(iMoment);
                end
                prob_aux(i_grp,ix1) = (1-mHat)* exp(mGridMoments*measureCoefficient')...
                    /normalization;
                
            end
            loglikes_hh(t,ix) = log(mean(prob_aux));
        end
    end
    
    loglikes2(i_beta) = loglikes2(i_beta)+sum(loglikes_hh(:));
    
    %% Save
    save(['loglikes_hh_beta' num2str(i_beta) '.mat'],'loglikes_hh')    
end

cd('../');