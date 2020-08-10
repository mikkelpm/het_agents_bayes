function likes = likelihood_micro(smooth_draw, t, data_micro, it, aaBar, mmu, ttau, mu_l, num_interp)

% Micro likelihood for household model

num_mom = 3;

likes = nan(length(data_micro,2),1);

for eepsilon=0:1 % For each employment status...
    
    ix = (data_micro(it,:,1)==eepsilon); % Households with given employment status

    % Collect distribution parameters
    moment = nan(num_mom,1);
    measureCoefficient = nan(num_mom,1);
    for im=1:num_mom
        moment(ii) = smooth_draw.(sprintf('%s%d%s%d', 'lag_moment_', eepsilon+1, '_', ii))(t);
        measureCoefficient(ii) = smooth_draw.(sprintf('%s%d%s%d', 'measureCoefficient_', eepsilon+1, '_', ii))(t);
    end
    mHat = smooth_draw.(sprintf('%s%d', 'mHat_', eepsilon+1))(t);
    
    % Compute normalization constant
    moment_aux = moment;
    moment_aux(1) = 0;
    g_log = @(a) measureCoefficient*((a-moment(1)).^((1:num_mom)')-moment_aux');
    lastwarn('');
    normalization = integral(@(a) exp(g_log(a)), aaBar, Inf);
    warnMsg = lastwarn;
    if ~isempty(warnMsg)
        disp(measureCoefficient);
        error('Improper asset density');
    end

    % Continuous part
    c = smooth_draw.w(t)*((1-eepsilon)*mmu+eepsilon*(1-ttau));
    R = smooth_draw.r(t);
    if R<=0
        warning('%s%8.4f', 'R=', R);
    end

    sigma2 = -2*mu_l;

    vals = linspace(min(log(data_micro(it,ix,2))),max(log(data_micro(it,ix,2))),num_interp); % Compute integral at these grid points for log income
    ints = zeros(1,num_interp);
    for i_in=1:num_interp
        ints(i_in) = integral(@(a) exp(g_log(a) ...
            -(0.5/sigma2)*((vals(i_in)-mu_l)-log(c+R*a)).^2 ...
            ), ...
            aaBar, Inf); % Compute integral at given grid point
    end
    likes(ix) = interp1(vals,ints,log(data_micro(it,ix,2)),'pchip'); % Cubic interpolation of integral between grid points
    likes(ix) = (likes(ix)./data_micro(it,ix,2)) * ((1-mHat)/(normalization*sqrt(2*pi*sigma2))); % Likelihood for continuous part

    % Point mass part
    likes(ix) = max(likes(ix) ...
        + (mHat/sqrt(2*pi*sigma2))*exp(-0.5/sigma2*(log(data_micro(it,ix,2))-log(c+R*aaBar)-mu_l).^2)./data_micro(it,ix,2), ...
        eps);
    
end

end