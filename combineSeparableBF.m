function A = combineSeparableBF(Ax, Ay)
% This function combines the evaluated basis function in Ax and Ay to spatial
% basis function under the assumption that they are separable.
% The lines in Ax and Ay must correspond to the same point, so that each rows of
% the combined matrix is defined as
% A(i,:) = kron(Ay(i,:), Ax(i,:));
%
% Input:
% Ax ... double [nxp], evaluated basis function at the x-coordinates
% Ay ... double [nxq], evaluated basis function at the y-coordinates
%
% Output:
% A ... double [nx(p*q)], combined matrix with evaluated, separable basis 
%                         functions

% Initialize matrix
A = zeros(size(Ax,1), size(Ax,2)*size(Ay,2));
% A= kron (Ay, Ax);

end