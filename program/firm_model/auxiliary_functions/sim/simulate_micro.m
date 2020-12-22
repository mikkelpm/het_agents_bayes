function simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, trunc_logn)

% Simulate household employment and income

global nMeasure ttheta nnu; % prodMin prodMax capitalMin capitalMax;

T_micro = length(ts_micro);
simul_data_micro = nan(T_micro,N_micro,2);

disp('Simulating micro data...');


for it=1:T_micro
    t = ts_micro(it);
    fprintf('t = %4d\n', t);
    
    % prepare the parameters
    nMeasure_all = (nMeasure+3)*nMeasure/2;
    moment = nan(1,nMeasure_all);
    for i_Measure = 1:nMeasure_all
        moment(i_Measure) = sim_struct.(['lag_moment_' num2str(i_Measure)])(t);
    end
    
    % draw (prod,logk)
    the_mean = moment(1:2);
    the_varcov = [moment(3) moment(4); moment(4) moment(5)];
    the_draws = mvnrnd(the_mean, the_varcov, N_micro);
    
    % transform to (logn, logk)
    c = log(nnu)+sim_struct.aggregateTFP(t)-sim_struct.logWage(t);
    simul_data_micro(it,:,:) = [(c+the_draws(:,1)+ttheta*the_draws(:,2))/(1-nnu) the_draws(:,2)];
    
end

% Truncate
simul_data_micro(repmat(simul_data_micro(:,:,1)<trunc_logn,1,1,2)) = nan;

end