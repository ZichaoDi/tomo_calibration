% fista_solver_instrumented.m - A FISTA implementation that logs performance metrics.

function [x, history] = fista_solver_instrumented(obj_f, obj_g, grad_f, prox_g, x0, Lf, params, max_iters, x_star)
% Purpose:
%   Solves min f(x) + g(x) using FISTA and returns a history of performance
%   metrics for analysis.
%
% Inputs:
%   obj_f:      Function handle for the smooth objective, e.g., @(x) 0.5*sum((A*x-b).^2).
%   obj_g:      Function handle for the non-smooth objective, e.g., @(x) lambda*norm(x,1).
%   grad_f:     Function handle for the gradient of f.
%   prox_g:     Function handle for the proximal operator of g.
%   x0:         The initial guess.
%   Lf:         The Lipschitz constant of grad_f.
%   lambda:     The regularization parameter for g.
%   max_iters:  The maximum number of iterations to run.
%   x_star:     (Optional) The ground truth solution for error calculation.
%   tol:        (Optional) Tolerance for early stopping.
%
% Outputs:
%   x:          The final solution.
%   history:    A struct containing vectors of performance metrics per iteration.

lambda = params.lambda;
tolerance = params.fine_tolerance;

%% 1. Handle Optional Inputs and Initialization
if nargin < 9
    x_star = []; % No ground truth provided
end

% Initialize iterates
x_k = x0;
y_k = x0;
t_k = 1;

% Pre-allocate history arrays for speed
history.iter      = zeros(max_iters, 1);
history.obj_val   = zeros(max_iters, 1);
history.grad_norm = zeros(max_iters, 1);
history.gt_error  = zeros(max_iters, 1);
history.x_change  = zeros(max_iters, 1);
history.y_k_hist  = zeros(numel(x0), max_iters); % <<< ADD THIS LINE


fprintf('--- Running Instrumented FISTA for %d iterations ---\n', max_iters);
tic; % START THE TIMER
%% 2. Main FISTA Loop
for k = 1:max_iters

    % Store old value for momentum update and convergence check
    x_k_old = x_k;

    
    
    
    % --- Core FISTA Algorithm Steps ---
    grad_y = grad_f(y_k);
    
    % Note: The lambda passed to prox_g must be scaled by the step size 1/Lf
    prox_step = prox_g(y_k - (1/Lf) * grad_y, lambda); %lambda / Lf

    x_k = max(0, prox_step);%non negativity of the elements
    
    t_k_old = t_k;
    t_k = (1 + sqrt(1 + 4*t_k_old^2)) / 2;
    
    y_k = x_k + ((t_k_old - 1) / t_k) * (x_k - x_k_old);

    history.iter(k) = k;
    history.y_k_hist(:, k) = y_k - x_k; % <<< ADD THIS LINE
    
    % 2) Objective function value: f(x_k) + g(x_k)
    history.obj_val(k) = obj_f(x_k) + obj_g(x_k);
    
    % 3) Norm of the Gradient Mapping (a better optimality measure for composite problems)
    grad_map = Lf * (y_k - x_k); % This is equivalent to G(y_k)
    history.grad_norm(k) = norm(grad_map);
    
    % 4) Error against ground truth (if provided)
    if ~isempty(x_star)
        history.gt_error(k) = norm(x_k - x_star);
    end
    
    % 5) Change from previous iterate
    history.x_change(k) = norm1(x_k - x_k_old)/numel(x_k);

    if history.x_change(k) < tolerance  && k > 1
        fprintf('Convergence reached at iteration %d.\n', k);
        break;
    end
    
    % --- Check for early convergence ---
    
    
end
endtime = toc;
%% 3. Trim History Arrays
% If the loop stopped early, trim the excess zeros from the history arrays.
if k < max_iters
    history.iter      = history.iter(1:k);
    history.obj_val   = history.obj_val(1:k);
    history.grad_norm = history.grad_norm(1:k);
    history.gt_error  = history.gt_error(1:k);
    history.x_change  = history.x_change(1:k);
    history.y_k_hist  = history.y_k_hist(:, 1:k); % <<< ADD THIS LINE
end

x = x_k; % Return the last iterate

history.elapsed_time = endtime; % STOP THE TIMER and store the result
% fprintf('--- FISTA Finished in %.4f seconds ---\n', history.elapsed_time);

end