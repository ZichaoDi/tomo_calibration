% calculate_smoothed_gradient.m - Calculates the gradient of the smoothed objective function

function grad_F_mu = calculate_smoothed_gradient(x_k, L_fine, s, lambda, params)
% Purpose:
%   Calculates the gradient of the smoothed version of the objective function, F_mu.
%   This is required for the decision logic (equation 3.4) in the main MAGMA loop.
%
%   The objective is F(x) = f(x) + g(x), where:
%     f(x) = (1/2) * ||L*x - s||^2  (smooth data fidelity term)
%     g(x) = lambda * ||x||_1        (non-smooth regularization term)
%
%   The smoothed version is F_mu(x) = f(x) + g_mu(x), where g_mu is a smooth
%   approximation of g(x), given by equation (2.3) in the paper.
%
%   The gradient is therefore: grad_F_mu = grad_f(x) + grad_g_mu(x).
%
% Inputs:
%   x_k:      The current iterate vector.
%   L_fine:   The fine-level system matrix (L or A).
%   s:        The sinogram (measurement) vector (b).
%   lambda:   The L1 regularization parameter.
%   params:   A struct containing algorithm parameters, including 'mu'.
%
% Output:
%   grad_F_mu: The gradient vector of the smoothed objective function.

% Unpack the smoothing parameter 'mu' from the params struct.
% 'mu' is a small positive constant that controls the degree of smoothing.
if isfield(params, 'mu')
    mu = params.mu;
else
    mu = 1e-6; % A reasonable default value
end

% Step 1: Calculate the gradient of the smooth part, grad_f(x).
% grad_f(x) = L' * (L*x - s)
grad_f = L_fine' * (L_fine * x_k - s);

% Step 2: Calculate the gradient of the smoothed non-smooth part, grad_g_mu(x).
% The formula for the gradient of the smoothed L1-norm is derived from
% equation (2.3) and given in the paper.
% It is applied element-wise.
% grad_g_mu(x)_j = (lambda * x_j) / sqrt(mu^2 + x_j^2)
% We use MATLAB's element-wise operators (.*, ./, .^) for efficiency.
grad_g_mu = (lambda * x_k) ./ sqrt(mu^2 + x_k.^2);

% Step 3: Combine the two gradients to get the final smoothed gradient.
grad_F_mu = grad_f + grad_g_mu;

end