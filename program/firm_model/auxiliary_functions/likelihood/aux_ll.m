function [the_loglike, the_loglike_macro, the_loglike_micro] = ...
         aux_ll(data_micro, ts_micro, ...
                num_smooth_draws, num_burnin_periods, ...
                trunc_logn, likelihood_type, ...
                M_, oo_, options_, ...
                update_param, varargin)
				
		% Different versions of full-information/macro likelihood, firm model
        
        global mat_suff;

        saveParameters;         % Save parameter values to files
        if update_param
            setDynareParameters;    % Update Dynare parameters in model struct
            compute_steady_state;   % Compute steady state, no need for parameters of agg dynamics
        end
        
        % Macro state variables used in micro likelihood
        smooth_vars = [{'logWage'; 'aggregateTFP'};
                       str_add_numbers('lag_moment_', 1:5)];
        
        % Parameters passed to micro likelihood function
        param = [nnu ttheta trunc_logn];
        
        % Log likelihood computation
        switch likelihood_type
            case 1 % Macro + full info micro
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, num_smooth_draws, ...
                                   M_, oo_, options_, ...
                                   data_micro, ts_micro, param, varargin{:});
            case 2 % Macro only
                [the_loglike, the_loglike_macro, the_loglike_micro] = ...
                    loglike_compute(strcat('simul', mat_suff, '.mat'), ...
                                   num_burnin_periods, smooth_vars, 0, ...
                                   M_, oo_, options_, ...
                                   [], [], param);
        end
        
end