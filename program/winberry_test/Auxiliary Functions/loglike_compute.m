function [loglike, loglike_macro, loglike_hh, smooth_draws] = loglike_compute(sim_macro, simul_data_hh, ts_hh, num_smooth_draws, num_burnin_periods, constr_tol, M_, oo_, options_)

% Compute likelihood for Krusell-Smith model

global aaBar nEpsilon nMeasure mmu ttau;


%% Macro likelihood and simulation smoother

% Determine variables to smooth
smooth_vars_aux1 = cell(nEpsilon,nMeasure);
smooth_vars_aux2 = cell(nEpsilon,nMeasure);
for i_Measure = 1:nMeasure
    for i_Epsilon = 1:nEpsilon
        smooth_vars_aux1{i_Epsilon,i_Measure} = ['moment_' num2str(i_Epsilon) '_' num2str(i_Measure)];
        smooth_vars_aux2{i_Epsilon,i_Measure} = ['measureCoefficient_' num2str(i_Epsilon) '_' num2str(i_Measure)];
    end
end
smooth_vars = char([{'w'; 'r'; 'mHat_1' ; 'mHat_2'}; smooth_vars_aux1(:); smooth_vars_aux2(:)]); %char(setxor(cellstr(M_.endo_names), options_.varobs)); % All endogenous variables except observed ones

% Run smoother and compute macro likelihood
[loglike_macro, ~, smooth_draws] = simulation_smoother(sim_macro, smooth_vars, num_smooth_draws, num_burnin_periods, M_, oo_, options_);


%% Household likelihood per period, given smoothing draws

T_hh = length(ts_hh);
loglikes_hh = nan(1,num_smooth_draws);

for i_draw = 1:num_smooth_draws
    
    the_loglikes_hh_draw = nan(1,T_hh);

    for it = 1:T_hh
    
        t = ts_hh(it);

        the_loglikes_hh_draw_t = nan(1,length(simul_data_hh(it,:,1)));

        for eepsilon = 0:1

            ix = simul_data_hh(it,:,1)==eepsilon;
            the_likes = nan(1,sum(ix));

            % assets
            asset_aux = (simul_data_hh(it,ix,2)-smooth_draws{i_draw}.w(t)*((1-eepsilon)*mmu+eepsilon*(1-ttau)))/(1+smooth_draws{i_draw}.r(t));

            % probability
            % prepare the parameters
            ix0 = abs(asset_aux-aaBar)<constr_tol;
            ix1 = ~ix0;
            mHat = smooth_draws{i_draw}.(['mHat_' num2str(eepsilon+1)])(t);
            moment = nan(1,nMeasure);
            measureCoefficient = nan(1,nMeasure);
            for i_Measure = 1:nMeasure
                moment(i_Measure) = smooth_draws{i_draw}.(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
                measureCoefficient(i_Measure) = smooth_draws{i_draw}.(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)])(t);
            end

            % constrained assets
            the_likes(ix0) = mHat;

            % unconstrained assets
            moment_aux = moment;
            moment_aux(1) = 0;
            g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux'));
            normalization = integral(g, aaBar, Inf);
            the_likes(ix1) = (1-mHat)*g(asset_aux(ix1))/normalization;

            the_loglikes_hh_draw_t(ix) = log(the_likes);

        end
        
        the_loglikes_hh_draw(it) = sum(the_loglikes_hh_draw_t);

    end
    
    loglikes_hh(i_draw) = sum(the_loglikes_hh_draw);
    
end


%% Sum log likelihood

log_max = max(loglikes_hh);
loglike_hh = log_max + log(mean(exp(loglikes_hh-log_max))); % Formula deals with underflow

loglike = loglike_macro + loglike_hh;

end