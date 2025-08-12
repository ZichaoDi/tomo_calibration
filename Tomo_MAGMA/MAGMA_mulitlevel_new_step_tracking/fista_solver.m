% fista_solver.m - A standard FISTA implementation.

function [x] = fista_solver(grad_f, prox_g, x0, Lf, lambda, max_iters)
% Purpose:
%   Solves min f(x) + g(x) using the Fast Iterative Shrinkage-Thresholding
%   Algorithm (FISTA).
%
% Inputs:
%   grad_f:     A function handle to compute the gradient of f, e.g., @(x) grad(x).
%   prox_g:     A function handle for the proximal operator of g, e.g., @(x,l) prox(x,l).
%   x0:         The initial guess.
%   Lf:         The Lipschitz constant of grad_f.
%   lambda:     The regularization parameter for g.
%   max_iters:  The number of iterations to run.
%
% Output:
%   x:          The approximate solution.

x_k = x0;
y_k = x0;
t_k = 1;

% tic;
for k = 1:max_iters
    
    % Store old value for momentum update
    x_k_old = x_k;
    
    % FISTA algorithm steps
    grad_y = grad_f(y_k);

    prox_step = prox_g(y_k - (1/Lf) * grad_y, lambda);
    
    % --- ENFORCE NON-NEGATIVITY ---
    % Project the result onto the non-negative set by setting all negative
    % elements to zero. This is the FISTA equivalent of your x(x<0)=0 step.
    x_k = max(0, prox_step);
    
    t_k_old = t_k;
    t_k = (1 + sqrt(1 + 4*t_k_old^2)) / 2;
    
    y_k = x_k + ((t_k_old - 1) / t_k) * (x_k - x_k_old);
    
end
% time = toc;
x = x_k; % Return the last iterate

end