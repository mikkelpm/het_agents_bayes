function mPoly = computeChebyshev(nPower,vGrid);

% Compute grid size
[nGrid,~] = size(vGrid);

% Create polynomial
mPoly = ones(nGrid,nPower);
mPoly(:,2) = vGrid;
for iPower = 3:nPower
	mPoly(:,iPower) = 2 * vGrid .* mPoly(:,iPower-1) - mPoly(:,iPower-2);
end