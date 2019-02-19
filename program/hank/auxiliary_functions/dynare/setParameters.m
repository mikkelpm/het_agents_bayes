% Sets parameter values

% Invariant idiosyncratic distribution
nz = length(vzGrid);
[zV,zD] = eig(mzTransition');
vzInvariant = zV(:,abs(diag(zD)-1)<tolerance);
assert(size(vzInvariant,2)==1);
vzInvariant = vzInvariant/sum(vzInvariant); % Invariant distribution
z_SS = vzGrid'*vzInvariant; % Steady state mean

% Share grid
nShare = length(vShareGrid);
assert(abs(vShareFraction'*vShareGrid-1)<1e-10); % Check that profit shares add up

% Approximation dimensions
nState = nz * nAssets;
nStateFine = nz * nAssetsFine;
nStateQuadrature = nz * nAssetsQuadrature;
nMeasureCoefficients = nz * nMeasure;

% Compute some steady state quantities, including for repr agent model
w_SS = (eepsilon-1)/eepsilon*A_SS;
r_RepSS = 1/bbeta-1;
N_RepSS = (((1-ttau)*w_SS/ppsi)/((1-r_RepSS*vvarthetaB+vvarthetaT)*A_SS-ttau*w_SS))^(1/(1+nnu))*z_SS;

% Bounds on grid space
assetsMin = bbBar;
assetsMax = 20 * (-vvarthetaB)*A_SS*N_RepSS; % Multiple of steady state borrowing in repr agent model
