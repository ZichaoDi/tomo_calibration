% prox_g.m - Proximal operator for the L1 norm

function x = prox_g(y, lambda_val)
% Purpose:
%   Implements the proximal operator for g(x) = lambda * ||x||_1.
%   This is the element-wise soft-shrinkage/thresholding operator.
%
%   Formula: x_i = sign(y_i) * max(|y_i| - lambda, 0) for each element i.
%
% Inputs:
%   y:          The input vector to be thresholded.
%   lambda_val: The thresholding parameter (lambda).
%
% Output:
%   x:          The thresholded vector.

% This single line performs the element-wise operation efficiently in MATLAB.
x = sign(y) .* max(abs(y) - lambda_val, 0);

end