clear all;


%% Initial Dynare run

dynare example1_new noclearall;         % Run Dynare once to process model file


%% Simulate data

% Settings
options_.periods = 100;                 % Simulate 100 periods
options_.drop = 0;                      % Don't drop any simulation periods

% Simulate
set_dynare_seed(20171215);              % Seed RNG
stoch_simul(var_list_);                 % Simulate

% Store simulated data
sim_vars = cellstr(options_.varobs)';   % Names of observed simulated
save('simul.mat', sim_vars{:});         % Save simulated data


%% Smoothing and likelihood

% Settings
options_.datafile = 'simul.mat';                % Data file
options_.periods = 0;                           % Don't simulate
options_.noprint = 1;                           % Don't print output after solving model
options_.mode_compute = 0;                      % Don't compute mode when "estimating"
options_.smoother = 1;                          % Return smoothed variables
options_.selected_variables_only = 1;           % Smooth selected variables (suppresses 
options_.smoothed_state_uncertainty = true;     % Return smoothing variance
options_.debug = true;                          % Returns likelihood value when running smoother
options_.filter_covariance = true;              % Returen [3D array]  (m*m*(T+1)) array of one-step ahead forecast error variance
%                                                 matrices (decision-rule order)

% Determine variables to smooth
smooth_vars = char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones

% Loop over alpha values
alphas = linspace(0.2, 0.5, 10); % alpha values to loop over
loglikes = zeros(size(alphas));

% tic
for i=1%:length(alphas) % For each alpha...

    fprintf('%s%4.2f\n', 'alpha=', alphas(i));
    set_param_value('alpha', alphas(i));                    % Set alpha
%     disp('stoch_simul')
%     stoch_simul(var_list_);                                 % Solve model with new parameters
%     disp('estimation')
%     disp(['ss ' num2str(oo_.steady_state')])
    dynare_estimation(smooth_vars);                         % Run smoother
    loglikes(i) = -oo_.likelihood_at_initial_parameters;    % Store log likelihood value

    %% draw smoothed state vars
    T = options_.nobs;
    
    v_state = oo_.dr.state_var;
    n_state = length(v_state);
    v_msure = dataset_info.missing.aindex{1}';
    n_msure = length(v_msure);
    
    T_mat = oo_.dr.ghx(oo_.dr.inv_order_var(oo_.dr.state_var)',:);
    F_t = oo_.Smoother.Variance(v_msure,v_msure,1:end-1);
    C_t = oo_.Smoother.Variance(v_state,v_msure,1:end-1);
    P_t_t1 = oo_.Smoother.Variance(v_state,v_state,1:end-1);
    
    % can be done in one step
    a_t_t = [];
    for i_state = 1:n_state
        a_t_t = [a_t_t eval(['oo_.UpdatedVariables.' M_.endo_names(oo_.dr.state_var(i_state),:)])];
    end
    a_t_t = a_t_t';
    J_t = nan(n_state,n_state,T);
    V_t = nan(n_state,n_state,T);
    for t = 1:T
        P_t_t = P_t_t1(:,:,t)-(C_t(:,:,t)/F_t(:,:,t))*C_t(:,:,t)';
        if t<T
            J_t(:,:,t) = P_t_t*T_mat'/P_t_t1(:,:,t+1);
            V_t(:,:,t) = P_t_t-J_t(:,:,t)*P_t_t1(:,:,t+1)*J_t(:,:,t)';
        else
            J_t(:,:,t) = 0;
            V_t(:,:,t) = P_t_t;
        end
    end
    
    N_draw = 200;
    state_draw = nan(n_state,T,N_draw);
    for i_draw = 1:N_draw % can change to parfor, can optimize through loops
        for t = T:-1:1
            if t<T
                mm = a_t_t(:,t)+J_t(:,:,t)*(state_draw(:,t+1,i_draw)-oo_.dr.ys(v_state)-T_mat*(a_t_t(:,t)-oo_.dr.ys(v_state)));
            else
                mm = a_t_t(:,t);
            end
            state_draw(:,t,i_draw) = mvnrnd(mm',(V_t(:,:,t)+V_t(:,:,t)')/2+1e-10*eye(n_state));
        end
    end
%%    
end
% disp(['ss ' num2str(oo_.steady_state')])
% t=toc