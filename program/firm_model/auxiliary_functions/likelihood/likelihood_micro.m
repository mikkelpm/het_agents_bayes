function likes = likelihood_micro(smooth_draw_t, data_micro_t, param)

    % Micro likelihood for firm model

    % Read parameters
    nnu = param(1);
    ttheta = param(2);
    trunc_logn = param(3);

    % Likelihood contribution
    c = log(nnu)+smooth_draw_t{1,'aggregateTFP'}-smooth_draw_t{1,'logWage'};
    jacob = 1-nnu;
    ix = ~isnan(data_micro_t(:,1));
    produ = jacob*data_micro_t(ix,1)-c-ttheta*data_micro_t(ix,2);

    the_mean = smooth_draw_t{1,{'lag_moment_1','lag_moment_2'}};
    the_varcov = [smooth_draw_t{1,{'lag_moment_3','lag_moment_4'}};
                  smooth_draw_t{1,{'lag_moment_4','lag_moment_5'}}];
    likes = jacob*mvnpdf([produ data_micro_t(ix,2)], the_mean, the_varcov);

    % Truncation term
    the_aux = [1 ttheta]';
    likes = likes/normcdf((the_mean*the_aux+c-(1-nnu)*trunc_logn)/sqrt(the_aux'*the_varcov*the_aux));

end