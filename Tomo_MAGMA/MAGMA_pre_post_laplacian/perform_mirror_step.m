% perform_mirror_step.m - Correctly handles fine (proximal) and coarse (gradient) steps.
function z_k_plus_1 = perform_mirror_step(z_k, x_k, alpha_k_plus_1, grad_at_xk, prox_op, lambda, is_fine_level)
% Purpose:
%   Performs the mirror descent update step (z_k -> z_{k+1}).
%   The logic depends on whether the level is fine (non-smooth) or coarse (smooth).
%
% Inputs:
%   z_k:            The current mirror point.
%   x_k:            The current query point.
%   alpha_k_plus_1: The step size for the mirror update.
%   grad_at_xk:     The gradient used for the update.
%                   - On fine level: This is grad(f(x_k)).
%                   - On coarse level: This is grad(F_H(x_k)).
%   prox_op:        The proximal operator.
%   lambda:         The regularization parameter.
%   is_fine_level:  A boolean (true/false) indicating the level type.

if is_fine_level
    % --- FINE LEVEL: Perform a true proximal step for mirror descent ---
    % This solves: min { 0.5*||z - z_k||^2 + <alpha*grad(f), z> + alpha*g(z) }
    % The solution is a proximal update on g.
    
    % The step-size for the prox is scaled by alpha.
    z_k_plus_1 = prox_op(z_k - alpha_k_plus_1 * grad_at_xk, alpha_k_plus_1 * lambda);
    
else
    % --- COARSE LEVEL: Perform a simple gradient descent step ---
    % The problem is fully smooth (g=0), so the mirror step simplifies.
    % Here, grad_at_xk is the gradient of the entire coarse objective F_H.
    
    z_k_plus_1 = z_k - alpha_k_plus_1 * grad_at_xk;
end

% Enforce non-negativity constraint, which is common to both cases.
z_k_plus_1 = max(0, z_k_plus_1);

end
