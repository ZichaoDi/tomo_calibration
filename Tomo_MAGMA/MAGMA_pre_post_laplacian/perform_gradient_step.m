% perform_gradient_step.m - Performs a single proximal gradient step

function y_k_plus_1 = perform_gradient_step(x_k, L_fine, s, lambda, Lf_fine)
% Purpose:
%   Performs a single proximal gradient step (also known as an ISTA/FISTA step).
%   This is used when a "gradient correction step" is chosen in the MAGMA algorithm.
%   The step consists of a standard gradient descent on the smooth part (f)
%   followed by the application of the proximal operator of the non-smooth part (g).
%
% Inputs:
%   x_k:      The current iterate vector.
%   L_fine:   The fine-level system matrix (L or A).
%   s:        The sinogram (measurement) vector (b).
%   lambda:   The L1 regularization parameter.
%   Lf_fine:  The Lipschitz constant of the gradient of f.
%
% Output:
%   y_k_plus_1: The result of the proximal gradient step.


% Step 1: Calculate the gradient of the smooth part of the objective.
% The smooth part is f(x) = (1/2) * ||L*x - s||^2
% The gradient is grad_f(x) = L' * (L*x - s)
grad_f = L_fine' * (L_fine * x_k - s);

% Step 2: Perform a standard gradient descent step.
% The step size is gamma = 1 / Lf_fine.
x_update = x_k - (1 / Lf_fine) * grad_f;

% Step 3: Apply the proximal operator for the L1 norm to the updated vector.
% The threshold parameter for the prox operator must also be scaled by the step size.
lambda_prox = lambda / Lf_fine;

% Call the soft-thresholding function we created earlier.
y_k_plus_1 = prox_g(x_update, lambda_prox);

end