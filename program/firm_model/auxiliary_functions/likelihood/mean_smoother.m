function oo_smooth = mean_smoother(data,xparam1,dataset_,dataset_info,M_,oo_,options_,bayestopt_,estim_params_)

    % Mean smoother

    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp,Trend,state_uncertainty,M_,oo_,options_,bayestopt_] = DsgeSmoother(xparam1,dataset_.nobs,data,dataset_info.missing.aindex,dataset_info.missing.state,M_,oo_,options_,bayestopt_,estim_params_);
    oo_smooth = store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty);
    
end