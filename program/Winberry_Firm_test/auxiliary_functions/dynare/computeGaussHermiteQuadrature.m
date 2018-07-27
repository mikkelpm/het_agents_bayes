function [grid,weight] = computeGaussHermiteQuadrature(order)

% This function computes the notes and weights for the Gauss - Hermite quadrature
% of order "order." 
%
% Inputs
%   (1) order: order of the quadrature
%
% Outputs
%   (1) grid: nodes to evaluate the function on
%   (2) weights: weight in the approximation
%
% Jung Sakong, February 10, 2016

% Compute polynomial recursively
temp_poly = zeros(order+1,order+1);
temp_poly(1,1) = 1;
temp_poly(2,2) = 2;
for ii = 3:(order+1)
	temp_poly(ii,1) = -temp_poly(ii-1,2);
	temp_poly(ii,2:end) = 2*temp_poly(ii-1,1:end-1) ...
		- ([2:1:order+1] .* [temp_poly(ii-1,3:end),0]);
end
the_poly = fliplr(temp_poly(end,:));
coef_before = fliplr(temp_poly(end-1,1:end-1)); % higher order coefficients first

% Solve for roots of the polynomial
grid = roots(the_poly);

% Compute weights
temp_powers = zeros(order,order);
for ii = 1:order
	temp_powers(:,ii) = grid.^(order-ii);
end
poly_before = repmat(coef_before,[order 1]) ...
	.* temp_powers;
weight = (2^(order-1)*factorial(order)*sqrt(pi))...
	./ ((order^2)*(sum(poly_before,2)).^2);