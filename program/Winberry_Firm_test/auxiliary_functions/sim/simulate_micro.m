function simul_data_micro = simulate_micro(sim_struct, ts_micro, N_micro, num_interp)

% Simulate household employment and income

global nMeasure ttheta nnu prodMin prodMax capitalMin capitalMax;

T_micro = length(ts_micro);
simul_data_micro = nan(T_micro,N_micro,2);

disp('Simulating micro data...');


for it=1:T_micro
    t = ts_micro(it);
    fprintf('t = %4d\n', t);
    
    % prepare the parameters
    nMeasure_all = (nMeasure+3)*nMeasure/2;
    moment = nan(1,nMeasure_all);
%     measureCoefficient = nan(1,nMeasure_all);
    for i_Measure = 1:nMeasure_all
        moment(i_Measure) = sim_struct.(['moment_' num2str(i_Measure)])(t-1);
%         measureCoefficient(i_Measure) = sim_struct.(['measureCoefficient_' num2str(i_Measure)])(t);
    end
    
    % draw (prod,logk)
    the_mean = moment(1:2);
    the_varcov = [moment(3) moment(4); moment(4) moment(5)];
    the_draws = mvnrnd(the_mean, the_varcov, N_micro);
    
    % transform to (logn, logk)
    c = log(nnu)+sim_struct.aggregateTFP(t)-sim_struct.logWage(t);
    simul_data_micro(it,:,:) = [(c+the_draws(:,1)+ttheta*the_draws(:,2))/(1-nnu) the_draws(:,2)];
    
    
%     % draw logy      
%     g = @(prod,logk) exp(g_kernel(prod,logk,moment,measureCoefficient,nMeasure));
%     normalization = integral2(g,-Inf,Inf,-Inf,Inf);
%     
%     ccdf_prod = @(p) integral2(@(prod,logk) g(prod,logk)/normalization, -Inf, p, -Inf, Inf);
%     
%     v_prod = nan(N_micro,1);
%     v_logk = nan(N_micro,1);
%     uu = rand(N_micro,2);
%     parfor i_micro = 1:N_micro
% %         if mod(i_micro,10) == 0
% %             fprintf('i_micro = %4d\n', i_micro);
% %         end
%         the_uu = uu(i_micro,:);
%         v_prod(i_micro) = fzero(@(p) ccdf_prod(p)-the_uu(1),moment(1));
%         the_norm = integral(@(logk) g(v_prod(i_micro),logk), -Inf, Inf);
%         the_ccdf_logk = @(lk) integral(@(logk) g(v_prod(i_micro),logk)/the_norm, -Inf, lk);
%         v_logk(i_micro) = fzero(@(lk) the_ccdf_logk(lk)-the_uu(2),moment(2));
%     end
%     
%     c = log(nnu)+sim_struct.aggregateTFP(t)-sim_struct.logWage(t);
%     simul_data_micro(it,:,:) = [(c+v_prod+ttheta*v_logk)/(1-nnu) v_logk];

    
%     the_vals = linspace(prodMin+ttheta*log(capitalMin),prodMax+ttheta*log(capitalMax),num_interp); % Compute integral at these grid points for log ouput
%     the_ints = zeros(1,num_interp);
%     for i_in=1:num_interp
%         the_ints(i_in) = integral2(@(x,prod) g(prod,(x-prod)/ttheta), ...
%             -Inf, the_vals(i_in), -Inf, Inf)/(ttheta*normalization);
%     end
%     pp = griddedInterpolant(the_vals,the_ints,'pchip'); % Cubic interpolation of integral between grid points
%     uu = rand(1,N_micro); % Random uniforms
%     uu_sort = sort(uu); % Sort uniforms to make it easier to do root finding
%     v_aux = nan(1,N_micro);
%     v_aux_prev = the_vals(1);
%     for i_micro = 1:N_micro
%         if uu_sort(i_micro)>the_ints(end)
%             v_aux(i_micro) = the_vals(end);
%         elseif uu_sort(i_micro)<the_ints(1)
%             v_aux(i_micro) = the_vals(1);
%         else
%             v_aux(i_micro) = fzero(@(x) pp(x)-uu_sort(i_micro), v_aux_prev); % Start root finding from previous solution
%         end
%         v_aux_prev = v_aux(i_micro); % Update previous solution
%     end
%     
%     % logy
%     simul_data_micro(it,:) = (v_aux+sim_struct.aggregateTFP(t)...
%         +nnu*(log(nnu)-log(sim_struct.wage(t))))/(1-nnu);
    
end
end