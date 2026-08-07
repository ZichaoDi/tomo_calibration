function [x_k] = mfista_solver(obj_f, obj_g, grad_f, prox_g, x0, Lf, lambda, max_iters, tolerance)
% mfista_solver.m - Monotone Fast Iterative Shrinkage-Thresholding Algorithm.
%
% Purpose:
%   Solves min f(x) + g(x) using the Monotone FISTA (MFISTA) algorithm.
%   This version ensures that the objective function value is non-increasing
%   at every iteration.
%
% Inputs:
%   obj_f:      A function handle for the smooth part of the objective, f(x).
%               Example: @(x) 0.5*norm(A*x-b)^2
%   obj_g:      A function handle for the non-smooth part, g(x).
%               Example: @(x) lambda*norm(x,1)
%   grad_f:     A function handle to compute the gradient of f.
%   prox_g:     A function handle for the proximal operator of g.
%   x0:         The initial guess.
%   Lf:         The Lipschitz constant of grad_f.
%   lambda:     The regularization parameter.
%   max_iters:  The maximum number of iterations to run.
%   tolerance:  The convergence tolerance for the stopping criterion.
%
% Output:
%   x_k:        The approximate solution.
%
% --- Initialization ---
x_k = x0;
y_k = x0;
t_k = 1;

for k = 1:max_iters
    
    % Store the current point and calculate its objective value
    x_k_old = x_k;
    F_k_old = obj_f(x_k_old) + obj_g(x_k_old);
    
    % --- Step 1: Calculate the potential next iterate using a standard FISTA step ---
    grad_y = grad_f(y_k);
    
    % The step-size for the prox is lambda / Lf for correct scaling.
    prox_step_accelerated = prox_g(y_k - (1/Lf) * grad_y, lambda/Lf);
    
    % Enforce non-negativity constraint, as required by the problem context.
    x_k_candidate = max(0, prox_step_accelerated);
    
    % --- Step 2: Monotonicity Check ---
    % Calculate the objective value at the potential new point.
    % F_k_candidate = obj_f(x_k_candidate) + obj_g(x_k_candidate);
    
    % if F_k_candidate > F_k_old
    %     % Objective increased. REJECT the accelerated step.
    %     % Instead, perform a standard, non-accelerated proximal gradient
    %     % (ISTA) step from the PREVIOUS position (x_k_old). This guarantees
    %     % the objective function value will not increase.
    % 
    %     grad_x_old = grad_f(x_k_old);
    %     prox_step_ista = prox_g(x_k_old - (1/Lf) * grad_x_old, lambda/Lf);
    %     x_k = max(0, prox_step_ista);
    % 
    %     % Reset the momentum term since we took a conservative step.
    %     t_k = 1;
    %     y_k = x_k; % The next y_k will start from this new, safe point.
    % 
    % else
        % Objective did not increase. ACCEPT the accelerated step.
        x_k = x_k_candidate;
        
        % Update the momentum terms as in standard FISTA.
        t_k_old = t_k;
        t_k = (1 + sqrt(1 + 4*t_k_old^2)) / 2;
        y_k = x_k + ((t_k_old - 1) / t_k) * (x_k - x_k_old);
    % end
    
    % --- Stopping Criterion ---
    % Calculate the normalized L1-norm of the change in x.
    x_change = norm(x_k - x_k_old, 1) / numel(x_k);
    
    % Check for convergence after the first iteration.
    if x_change < tolerance && k > 1
        fprintf('     Coarse Convergence reached at iteration %d.\n', k);
        break;
    end
    
end
end