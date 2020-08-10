function likes = likelihood_micro(smooth_draw, t, data_micro, it, nnu, ttheta, trunc_logn)

    % Micro likelihood for firm model
    
    c = log(nnu)+smooth_draw.aggregateTFP(t)-smooth_draw.logWage(t);
    jacob = 1-nnu;
    produ = jacob*data_micro(it,:,1)-c-ttheta*data_micro(it,:,2);

    the_mean = [smooth_draw.lag_moment_1(t) smooth_draw.lag_moment_2(t)];
    the_varcov = [smooth_draw.lag_moment_3(t) smooth_draw.lag_moment_4(t);
                  smooth_draw.lag_moment_4(t) smooth_draw.lag_moment_5(t)];
    likes = jacob*mvnpdf([produ' data_micro(it,:,2)'], the_mean, the_varcov);

    % Truncation term
    the_aux = [1 ttheta]';
    likes = likes/normcdf((the_mean*the_aux+c-(1-nnu)*trunc_logn)/sqrt(the_aux'*the_varcov*the_aux));

end