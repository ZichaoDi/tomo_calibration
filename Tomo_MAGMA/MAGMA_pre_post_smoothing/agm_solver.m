function [x, history] = agm_solver(obj_f, obj_g, grad_f, prox_op, Lf_current, s, x0, params, x_star)
%AGM_SOLVER An accelerated gradient method 
%of the user's MAGMA solver, with the coarse-step logic removed.

    % --- Unpack Parameters from params struct ---
    lambda = params.lambda;
    T = params.T;
    fine_tolerance = params.fine_tolerance;
    
    % --- Initialization (exactly like magma_solver) ---
    x_k = x0;
    y_k = x0; 
    z_k = x0;
    eta_k = Lf_current;
    alpha_k = 0; 
    
    % --- History Allocation ---
    history.iter = zeros(T, 1);
    history.obj_val = zeros(T, 1);
    history.gt_error = zeros(T, 1);
    history.x_change = zeros(T, 1);
    % Added history for plotting the vectors, as requested
    history.y_k_hist = zeros(numel(x0), T);

    fprintf('Starting AGM Solver (Strict MAGMA Structure, T=%d)...\n', T);

    % --- Main AGM Loop ---
    for k = 1:T
        x_k_old = x_k;
        
        % --- Core AGM Steps (Kept exactly from magma_solver) ---
        if alpha_k == 0, t_k = 1; else, t_k = 1 / (alpha_k * eta_k); end
        t_k = min(max(t_k, 0), 1);
        x_k = (1 - t_k) * y_k + t_k * z_k;
        
        % --- This is the Gradient Step logic from magma_solver's 'else' branch ---
        grad_at_xk = grad_f(x_k);
        y_k_plus_1 = prox_op(x_k - (1/Lf_current) * grad_at_xk, lambda/Lf_current);
        y_k_plus_1 = max(0, y_k_plus_1);
        
        % --- Parameter updates (Kept exactly from magma_solver) ---
        eta_k_plus_1 = Lf_current;
        alpha_k_plus_1 = (k + 2) / (2 * Lf_current);
        
        % --- Mirror Descent (Kept exactly from magma_solver) ---
        
        z_k_plus_1 = perform_mirror_step(z_k, x_k, alpha_k_plus_1, grad_at_xk, prox_op, lambda, true);
        
        history.y_k_hist(:, k) = y_k_plus_1 - y_k; % Save the y_k iterate for plotting

        % --- State Updates ---
        y_k = y_k_plus_1; 
        z_k = z_k_plus_1;
        eta_k = eta_k_plus_1; 
        alpha_k = alpha_k_plus_1;
        
        % --- Record History (Kept exactly as in magma_solver) ---
        history.iter(k) = k;
        history.obj_val(k) = obj_f(x_k) + obj_g(x_k); % Logging based on x_k
        history.gt_error(k) = norm(x_k - x_star);    % Logging based on x_k
        history.x_change(k) = norm1(x_k - x_k_old) / numel(x_k); % Logging based on x_k
        

        % --- Convergence Check (Kept exactly as in magma_solver) ---
        if history.x_change(k) < fine_tolerance && k > 1
            fprintf('Convergence reached at AGM iteration %d.\n', k);
            % Trim history arrays
            history.iter = history.iter(1:k);
            history.obj_val = history.obj_val(1:k);
            history.gt_error = history.gt_error(1:k);
            history.x_change = history.x_change(1:k);
            history.y_k_hist = history.y_k_hist(:, 1:k);
            break;
        end
    end
    
    x = y_k; % Final result is y_k, same as in magma_solver
end