function likes = likelihood_micro(smooth_draw_t, data_micro_t, param)

% Micro likelihood for household model

% Read parameters
aaBar = param(1);
mmu = param(2);
ttau = param(3);
mu_l = param(4);
num_mom = param(5);
num_interp = param(6);

% Likelihood contribution

likes = nan(size(data_micro_t,1),1);

for eepsilon=0:1 % For each employment status...
    
    ix = (data_micro_t(:,1)==eepsilon); % Households with given employment status

    % Collect distribution parameters
    moment = smooth_draw_t{1,str_add_numbers(sprintf('%s%d%s', 'lag_moment_', eepsilon+1, '_'), 1:num_mom)};
    measureCoefficient = smooth_draw_t{1,str_add_numbers(sprintf('%s%d%s', 'measureCoefficient_', eepsilon+1, '_'), 1:num_mom)};
    mHat = smooth_draw_t{1,sprintf('%s%d', 'lag_mHat_', eepsilon+1)};
    
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
    c = smooth_draw_t{1,'w'}*((1-eepsilon)*mmu+eepsilon*(1-ttau));
    R = smooth_draw_t{1,'r'};
    if R<=0
        warning('%s%8.4f', 'R=', R);
    end

    sigma2 = -2*mu_l;

    vals = linspace(min(log(data_micro_t(ix,2))),max(log(data_micro_t(ix,2))),num_interp); % Compute integral at these grid points for log income
    ints = zeros(1,num_interp);
    for i_in=1:num_interp
        ints(i_in) = integral(@(a) exp(g_log(a) ...
            -(0.5/sigma2)*((vals(i_in)-mu_l)-log(c+R*a)).^2 ...
            ), ...
            aaBar, Inf); % Compute integral at given grid point
    end
    likes(ix) = interp1(vals,ints,log(data_micro_t(ix,2)),'pchip'); % Cubic interpolation of integral between grid points
    likes(ix) = (likes(ix)./data_micro_t(ix,2)) * ((1-mHat)/(normalization*sqrt(2*pi*sigma2))); % Likelihood for continuous part

    % Point mass part
    likes(ix) = max(likes(ix) ...
        + (mHat/sqrt(2*pi*sigma2))*exp(-0.5/sigma2*(log(data_micro_t(ix,2))-log(c+R*aaBar)-mu_l).^2)./data_micro_t(ix,2), ...
        eps);
    
end

end