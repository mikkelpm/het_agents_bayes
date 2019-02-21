function sim_data_micro = simulate_micro(sim_struct, ts_micro, N_micro)

% Simulate household employment and income

load('economicParameters', 'bbBar', 'vzInvariant', 'vShareGrid', 'vShareFraction');
load('approximationParameters', 'nz', 'nMeasure', 'nShare');

T_micro = length(ts_micro);
sim_data_micro = nan(T_micro,N_micro,2);

disp('Simulating micro data...');

% Idiosyncratic productivity indexes
z_draws = reshape(sum(rand(T_micro*N_micro,1)<cumsum(vzInvariant)',2),T_micro,N_micro);

% Household type
sim_data_micro(:,:,1) = reshape(vShareGrid(sum(rand(T_micro*N_micro,1)<cumsum(vShareFraction)',2)),T_micro,N_micro);

% Assets
for i_t=1:T_micro % For each period with micro data...
    
    t = ts_micro(i_t);
    
    fprintf('t = %4d\n', t);
    
    for i_z = 1:nz % For each productivity level...
        
        for i_Share = 1:nShare % For each household type...
            
            ix = find((z_draws(i_t,:)==i_z) & (sim_data_micro(i_t,:,1)==vShareGrid(i_Share))); % Households with this productivity and of this type

            % Prepare the parameters
            mHat = sim_struct.(sprintf('mHat_%d_%d',i_z,i_Share))(t-1);
            moment = nan(1,nMeasure);
            measureCoefficient = nan(1,nMeasure);
            for i_Measure = 1:nMeasure
                moment(i_Measure) = sim_struct.(sprintf('moment_%d_%d_%d',i_z,i_Measure,i_Share))(t-1);
                measureCoefficient(i_Measure) = sim_struct.(sprintf('measureCoefficient_%d_%d_%d',i_z,i_Measure,i_Share))(t); % Different time subscript - cf. equations_polynomials.mod
            end

            % Identify constrained households
            n_ix = length(ix);
            ix0 = rand(1,n_ix)<mHat; % Constrained households
            ix1 = ~ix0; % Unconstrained households

            % Set constrained assets
            sim_data_micro(i_t,ix(ix0),2) = bbBar;

            % Draw unconstrained assets
            moment_aux = moment;
            moment_aux(1) = 0;
            g = @(b) exp(measureCoefficient*((b-moment(1)).^((1:nMeasure)')-moment_aux'));
            normalization = integral(g, bbBar, Inf);

            n_ix1 = sum(ix1);
            bb = nan(1,n_ix1);
            uu = rand(1,n_ix1);
            ccdf = @(b) integral(@(x) g(x)/normalization,bbBar,b);
            for i_ix1 = 1:n_ix1
                bb(i_ix1) = fzero(@(b) ccdf(b)-uu(i_ix1),moment(1));
            end
            sim_data_micro(i_t,ix(ix1),2) = bb;
            
        end
        
    end
    
end

end