function simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro)

% Simulate household employment and income

global nMeasure ttheta nnu;

T_micro = length(ts_micro);
simul_data_micro = nan(T_micro,N_micro);

disp('Simulating micro data...');


for it=1:T_micro
    t = ts_micro(it);
    fprintf('t = %4d\n', t);
    
    % prepare the parameters
    nMeasure_all = (nMeasure+3)*nMeasure/2;
    moment = nan(1,nMeasure_all);
    measureCoefficient = nan(1,nMeasure_all);
    for i_Measure = 1:nMeasure_all
        moment(i_Measure) = sim_struct.(['moment_' num2str(i_Measure)])(t-1);
        measureCoefficient(i_Measure) = sim_struct.(['measureCoefficient_' num2str(i_Measure)])(t-1);
    end
    
    % draw prod and logk      
    g = @(prod,logk) exp(g_kernel(prod,logk,moment,measureCoefficient));
    normalization = integral2(g,-Inf,Inf,-Inf,Inf);
    
    ccdf_prod = @(p) integral2(@(prod,logk) g(prod,logk)/normalization, -Inf, p -Inf, Inf);
    
    v_prod = nan(1,N_micro);
    v_logk = nan(1,N_micro);
    uu = rand(2,N_micro);
    for i_micro = 1:n_micro
        v_prod(i_micro) = fzero(@(p) ccdf_prod(p)-uu(1,i_micro),moment(1));
        ccdf_logk = @(lk) integral(@(logk) g(v_prod(i_micro),logk)/normalization, -Inf, lk) ...
            /integral(@(logk) g(v_prod(i_micro),logk)/normalization, -Inf, Inf);
        v_logk(i_micro) = fzero(@(lk) ccdf_logk(lk)-uu(2,i_micro),moment(2));
    end
    
    % logy
    simul_data_micro(it,:) = (v_prod+sim_struct.aggregateTFP(t)+ttheta*v_logk...
        +nnu*(log(nnu)-log(sim_struct.wage(t))))/(1-nnu);
    
end
end