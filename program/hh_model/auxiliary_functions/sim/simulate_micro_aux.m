function simul_data_micro_aux = simulate_micro_aux(sim_struct, ts_micro, N_micro)

% Simulate household employment and normalized income

global nMeasure aggEmployment aaBar mmu ttau;

T_micro = length(ts_micro);
simul_data_micro_aux = nan(T_micro,N_micro,3);

disp('Simulating household data...');

% Employment outcomes
simul_data_micro_aux(:,:,1) = (rand(T_micro,N_micro)<aggEmployment)+0;

% Income outcomes
for it=1:T_micro
    t = ts_micro(it);
    fprintf('t = %4d\n', t);
    for eepsilon = 0:1
        ix = simul_data_micro_aux(it,:,1)==eepsilon;

        % Prepare the parameters
        mHat = sim_struct.(['mHat_' num2str(eepsilon+1)])(t-1);
        moment = nan(1,nMeasure);
        measureCoefficient = nan(1,nMeasure);
        for i_Measure = 1:nMeasure
            moment(i_Measure) = sim_struct.(['lag_moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
            measureCoefficient(i_Measure) = sim_struct.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
        end

        % Draw asset
        n_ix = sum(ix);
        asset_aux = nan(1,n_ix);
        ix0 = rand(1,n_ix)<mHat;
        ix1 = ~ix0;

        % Constrained assets
        asset_aux(ix0) = aaBar;

        % Unconstrained assets
        moment_aux = moment;
        moment_aux(1) = 0;
        g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux'));
        normalization = integral(g, aaBar, Inf);

        n_ix1 = sum(ix1);
        aa = nan(1,n_ix1);
        uu = rand(1,n_ix1);
        ccdf = @(a) integral(@(x) g(x)/normalization,aaBar,a);
        for i_ix1 = 1:n_ix1
            aa(i_ix1) = fzero(@(a) ccdf(a)-uu(i_ix1),moment(1));
        end
        asset_aux(ix1) = aa;

        R = sim_struct.r(t);
        simul_data_micro_aux(it,ix,2) = sim_struct.w(t)*((1-eepsilon)*mmu+eepsilon*(1-ttau))+R*asset_aux;
        simul_data_micro_aux(it,ix,3) = asset_aux;
    end
end

end