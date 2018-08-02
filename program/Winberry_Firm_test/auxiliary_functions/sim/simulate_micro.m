function simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, num_interp)

% Simulate household employment and income

global nMeasure ttheta nnu prodMin prodMax capitalMin capitalMax;

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
    
    % draw logy      
    g = @(prod,logk) exp(g_kernel(prod,logk,moment,measureCoefficient,nMeasure));
    normalization = integral2(g,-Inf,Inf,-Inf,Inf);
    
    the_vals = linspace(prodMin+ttheta*log(capitalMin),prodMax+ttheta*log(capitalMax),num_interp); % Compute integral at these grid points for log ouput
    the_ints = zeros(1,num_interp);
    for i_in=1:num_interp
        the_ints(i_in) = integral(@(prod) g(prod,(the_vals(i_in)-prod)/ttheta)/normalization, ...
            -Inf, Inf);
    end
    pp = pchip(the_vals,the_ints); % Cubic interpolation of integral between grid points
    uu = rand(1,N_micro);  
    v_aux = nan(1,N_micro);  
    for i_micro = 1:N_micro
        ix = find(uu(i_micro)<the_ints,1,'first');
        if isempty(ix)
            v_aux(i_micro) = the_vals(end);
        elseif ix == 1
            v_aux(i_micro) = the_vals(1);
        else
            ix_1 = ix-1;
            v_aux(i_micro) = fzero(@(x) pp.coefs(ix_1,:)*((x-the_vals(ix_1)).^(3:-1:0)')-uu(i_micro),...
                (the_vals(ix_1)+the_vals(ix))/2);
        end
    end
    
    % logy
    simul_data_micro(it,:) = (v_aux+sim_struct.aggregateTFP(t)...
        +nnu*(log(nnu)-log(sim_struct.wage(t))))/(1-nnu);
    
end
end