function A = evalPolynomial(x, p)
% This function evaluates the monomial basis of a polynomial of degree p at the 
% points x, so that
%
% A = [x.^0 x.^1 x.^2 ... x.^p].
% 
% Input:
% x ... double [nx1], evaluation points
% p ... double, degree of the polynomial
%
% Output:
% A ... double [nx(p+1)], matrix with evaluated basis function at the points x.
% Initialize matrix

% A = zeros(length(x), p+1);
A = bsxfun(@power, x, 0:p);



end