n_ix = 1000;
asset_aux = nan(2,n_ix);
for eepsilon = 0:1
% prepare the parameters

moment = nan(1,nMeasure);
measureCoefficient = nan(1,nMeasure);
for i_Measure = 1:nMeasure
    moment(i_Measure) = eval(['moment_' num2str(eepsilon+1) '_' num2str(i_Measure)]);
    measureCoefficient(i_Measure) = eval(['measureCoefficient_' num2str(eepsilon+1) '_' num2str(i_Measure)]);
end

% draw asset

% unconstrained assets
moment_aux = moment;
moment_aux(1) = 0;
g = @(a) exp(measureCoefficient*((a-moment(1)).^((1:nMeasure)')-moment_aux'));
normalization = integral(g, aaBar, Inf);

aa = nan(1,n_ix);
uu = rand(1,n_ix);
ccdf = @(a) integral(@(x) g(x)/normalization,aaBar,a);
for i_ix = 1:n_ix
    aa(i_ix) = fzero(@(a) ccdf(a)-uu(i_ix),moment(1));
end
asset_aux(eepsilon+1,:) = aa;
end