%clear all;


%% Settings

burnin = 2e3;
acf_lags = 200;


%% MCMC diagnostics

load ./dynamicModel/metropolis/dynamicModel_mh1_blck1
post_draws = x2;
params = {'\rho_{\epsilon}', '\sigma_{\epsilon}'};

% Trace plot
figure('Units', 'normalize', 'Position', [0.2 0.1 0.6 0.8]);
for j=1:2
    subplot(2,1,j);
    plot(post_draws(:,j), '-k');
    the_ylim = ylim;
    line(burnin*[1 1], the_ylim, 'LineStyle', '--', 'Color', 'k');
    ylim(the_ylim);
    title(params{j}, 'FontSize', 14);
end

% ACF
figure('Units', 'normalize', 'Position', [0.2 0.1 0.6 0.8]);
for j=1:2
    subplot(2,1,j);
    the_acf = autocorr(post_draws(burnin+1:end,j), acf_lags);
    plot(0:acf_lags, the_acf, '-k', 'LineWidth', 2);
    ylabel('');
    the_xlim = xlim;
    line(the_xlim, [0 0], 'Color', 'k');
    xlim(the_xlim);
    title(params{j}, 'FontSize', 14);
end


%% Posterior densities

params_truth = [0.53 0.0364];

figure('Units', 'normalize', 'Position', [0.1 0.2 0.8 0.6]);
for j=1:2
    subplot(1,2,j);
    [the_f,the_xi] = ksdensity(post_draws(burnin+1:end,j));
    plot(the_xi, the_f, '-k', 'LineWidth', 2);
    the_ylim = ylim;
    line(params_truth(j)*[1 1], the_ylim, 'LineStyle', '--', 'Color', 'k');
    ylim(the_ylim);
    title(params{j}, 'FontSize', 14);
end

