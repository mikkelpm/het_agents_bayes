% Sets parameter values

% Invariant idiosyncratic distribution
nz = length(vzGrid);
[zV,zD] = eig(mzTransition);
vzInvariant = zV(abs(diag(zD)-1)<tolerance); % Invariant distribution
assert(size(vzInvariant,2)==1);
zzss = vzGrid'*vzInvariant; % Steady state mean

% Approximation dimensions
nState = nEpsilon * nAssets;
nStateFine = nz * nAssetsFine;
nStateQuadrature = nz * nAssetsQuadrature;
nMeasureCoefficients = nz * nMeasure;

% Bounds on grid space
assetsMin = aaBar;
assetsMax = 3 * kRepSS;
