function simul_data_hh = simulate_hh(sim_struct, ts_hh, N_hh)

% Simulate household employment and income

global nMeasure aggEmployment aaBar mmu ttau;

T_hh = length(ts_hh);
simul_data_hh = nan(T_hh,N_hh,2);

disp('Simulating household data...');

% Employment outcomes
simul_data_hh(:,:,1) = (rand(T_hh,N_hh)<aggEmployment)+0;

% Income outcomes
for it=1:T_hh
    t = ts_hh(it);
    fprintf('t = %4d\n', t);
    for eepsilon = 0:1
        ix = simul_data_hh(it,:,1)==eepsilon;

        % prepare the parameters
        mHat = sim_struct.(['mHat_' num2str(eepsilon+1)])(t-1);
        moment = nan(1,nMeasure);
        measureCoefficient = nan(1,nMeasure);
        for i_Measure = 1:nMeasure
            moment(i_Measure) = sim_struct.(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t-1);
            measureCoefficient(i_Measure) = sim_struct.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
        end

        % draw asset
        n_ix = sum(ix);
        asset_aux = nan(1,n_ix);
        ix0 = rand(1,n_ix)<mHat;
        ix1 = ~ix0;

        % constrained assets
        asset_aux(ix0) = aaBar;

        % unconstrained assets
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
        simul_data_hh(it,ix,2) = sim_struct.w(t)*((1-eepsilon)*mmu+eepsilon*(1-ttau))+R*asset_aux;
    end
end

end