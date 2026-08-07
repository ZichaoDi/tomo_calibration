% estimate_lipschitz_power.m - UPDATED to use a pre-calculated L'*L matrix

function Lf_estimate = estimate_lipschitz_power(L)
% Purpose:
%   Instantly estimates an upper bound for the Lipschitz constant (largest 
%   eigenvalue of LtL) using the matrix infinity-norm. This is derived
%   from the Gershgorin Circle Theorem.
%
%   This method is extremely fast but may overestimate the true value.
%
% Inputs:
%   LtL: The pre-calculated matrix (L' * L).
%
% Output:
%   Lf_estimate: An upper bound on the largest eigenvalue of LtL.

% The infinity norm ('inf') for a symmetric matrix computes the maximum 
% absolute row sum, which is a guaranteed upper bound on the max eigenvalue.
Lf_estimate = norm(L, 1) * norm(L, 'inf');

end