function [xEst, lEst] = linearGMM(A, l)
% This function solves a linear Gauss-Markov model under the assumption of 
% independent and identically distributed observations.
% 
% Input:
% A ... [nxm] double, design matrix
% l ... [nx1] double, vector of observations
 
% Output:
% xEst ... [mx1] double, estimated vector of unknowns
% lEst ... [nx1] double, vector of adjusted observation
% 

% set options for the solver of linsolve
 opts.SYM = true;
 opts.POSDEF = true;

% solve normal equations
xEst = linsolve( A'*A,A'*l ,opts);

% adjusted observation
lEst = A*xEst;



end

