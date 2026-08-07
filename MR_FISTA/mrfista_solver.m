function [x_k, history] = mrfista_solver(obj_f_k, obj_g_k, grad_f_k, prox_op_k, ...
                                        x0, params, level_no, x_star, v_parent)
% mrfista_solver: Implements the Multilevel FISTA algorithm recursively.
%
% This function is based on the logic from SolveMRFISTA.m, but adapted
% to the user's framework. It combines FISTA updates with a multilevel
% coarse correction, governed by conditions from the MPGM paper.

    % --- 1. Unpack Parameters and Operators for Current Level ---
    T = params.T;
    if level_no > 1
        T = floor(T * 1.5); % More iterations on coarser levels
    end
    
    lambda = params.lambda; mu = params.mu;
    kappa = params.kappa; theta = params.theta; Kd = params.Kd;
    s = params.s;
    is_fine_level = (level_no == 1);
    
    % Get operators for the *current* level
    L_k = params.L_level_full{level_no};
    Lf_k = params.Lf_level_full{level_no}; % Lipschitz constant
    
    num_levels = length(params.L_level_full);
    
    if level_no < num_levels
        R_k = params.R_level_full{level_no}; % Restriction R: k -> k+1
        P_k = params.P_level_full{level_no}; % Prolongation P: k+1 -> k
        L_coarse = params.L_level_full{level_no+1};
    end
    
    % --- 2. Initialization ---
    x_k = x0;
    x_k_old = x0;
    y_k = x0;
    t_k = 1;
    t_k_old = 1;
    
    q = 0; % Consecutive gradient step counter
    x_tilde = x0; % Last point where coarse correction was done
    
    % --- 3. History & Work Tracking ---
    history = [];
    total_equivalent_work = 0;
    my_level_factor = params.N_vec(level_no) / params.N_vec(1);
    
    if is_fine_level
        history.iter = zeros(T, 1);
        history.obj_val = zeros(T, 1);
        history.gt_error = zeros(T, 1);
        history.cumulative_work = zeros(T, 1);
        fprintf('Starting MR-FISTA on Fine Level 1 (T=%d)...\n', T);
    end

    % --- 4. Main MR-FISTA Loop ---
    for k = 1:T
        
        % --- FISTA Extrapolation Step ---
        y_k = x_k + (t_k_old - 1)/t_k * (x_k - x_k_old);
        
        % --- Calculate Gradients at y_k ---
        grad_f_y = grad_f_k(y_k);
        % Smoothed gradient of g(x) = lambda * ||x||_1
        grad_g_mu_y = (lambda * y_k) ./ sqrt(mu^2 + y_k.^2);
        % Gradient of the smoothed objective F_mu(y_k)
        grad_F_mu_y = grad_f_y + grad_g_mu_y;

        do_coarse_step = false;
        if level_no < num_levels
            % --- Decision Block (Check Coarse Correction Conditions) ---
            % Condition 1: Gradient is not "too small" on the coarse grid
            grad_F_mu_y_restricted = R_k * grad_F_mu_y;
            cond1 = norm(grad_F_mu_y_restricted) > kappa * norm(grad_F_mu_y);
            
            % Condition 2: Not too close to the last correction point
            % (using relative distance as in SolveMRFISTA)
            cond2 = (norm(y_k - x_tilde) > theta * norm(x_tilde)) || (q < Kd);
            
            do_coarse_step = cond1 && cond2;
        end
        
        % --- 5. Perform Step (Coarse Correction or Gradient) ---
        if do_coarse_step
            % --- 5a. Coarse Correction Step ---
            if is_fine_level, fprintf('  k=%d: Coarse Step\n', k); end
            
            % --- Define Coarse Problem ---
            x_H0 = R_k * y_k; % Coarse initial guess
            
            % Coarse gradient of f(x) at x_H0
            grad_f_coarse_at_H0 = L_coarse' * (L_coarse * x_H0 - s);
            % Coarse gradient of smoothed g(x) at x_H0
            grad_g_mu_coarse_at_H0 = (lambda * x_H0) ./ sqrt(mu^2 + x_H0.^2);
            % Full coarse smoothed gradient at x_H0
            grad_F_mu_coarse_at_H0 = grad_f_coarse_at_H0 + grad_g_mu_coarse_at_H0;
            
            % --- Calculate Coherency Term (v_H) ---
            % This is v = grad_0 - grad_lp1 from SolveMRFISTA
            v_H_local = R_k * grad_F_mu_y - grad_F_mu_coarse_at_H0;
            v_H_total = v_H_local + R_k * v_parent; 
            
            % --- Define Coarse Objective Functions (with v_H) ---
            obj_f_coarse  = @(x_H) 0.5 * norm(L_coarse * x_H - s)^2 + dot(v_H_total, x_H);
            obj_g_coarse  = @(x_H) lambda * norm(x_H, 1);
            grad_f_coarse = @(x_H) L_coarse' * (L_coarse * x_H - s) + v_H_total;
            prox_op_coarse = @prox_g;
            
            % --- Recursive Call to Coarse Solver ---
            coarse_params = params;
            coarse_params.T = params.NH; % Use NH iterations
            
            [x_H_star, history_coarse] = mrfista_solver( ...
                obj_f_coarse, obj_g_coarse, grad_f_coarse, prox_op_coarse, ...
                x_H0, coarse_params, level_no + 1, [], v_H_total ...
            );
        
            % --- Calculate Coarse Correction Direction ---
            % d_k = P * (x_H0 - x_H_star), as in SolveMRFISTA
            d_k = P_k * (x_H0 - x_H_star);
            
            % --- Line Search (Backtracking) ---
            % Find step 's_k' such that F(y_k - s_k*d_k) satisfies Armijo
            % We use the *smoothed* objective for the line search, per paper
            s_k = 1.0;
            F_mu_y = calculate_smoothed_objective(y_k, obj_f_k, lambda, mu);
            descent_term = dot(grad_F_mu_y, d_k);
            
            if descent_term < 0 % Not a descent direction
                fprintf('  Warning: Coarse direction not a descent direction. Skipping step.\n');
                x_k_plus_1 = y_k; % Do nothing
                work_this_iter = my_level_factor * 1; % Failed step still costs 1
            else
                while s_k > 1e-6
                    x_candidate = y_k - s_k * d_k;
                    F_mu_candidate = calculate_smoothed_objective(x_candidate, obj_f_k, lambda, mu);
                    
                    if F_mu_candidate <= F_mu_y - params.c * s_k * descent_term
                        break; % Armijo condition satisfied
                    end
                    s_k = s_k * params.tau;
                end
                x_k_plus_1 = x_candidate;
                work_this_iter = history_coarse.cumulative_work(end) + my_level_factor * 1;
            end
            
            q = 0; % Reset gradient step counter
            x_tilde = y_k; % Store this point
            
        else
            % --- 5b. Standard FISTA Step (Gradient Step) ---
            if is_fine_level && mod(k,10) == 0, fprintf('  k=%d: Grad Step\n', k); end
            
            % Standard FISTA proximal gradient step from point y_k
            x_k_plus_1 = prox_op_k(y_k - (1/Lf_k) * grad_f_y, lambda/Lf_k);
            
            q = q + 1; % Increment gradient step counter
            work_this_iter = my_level_factor * 1; % 1 unit of work
        end
        
        % --- 6. FISTA State Updates ---
        t_k_plus_1 = 0.5 * (1 + sqrt(1 + 4*t_k^2));
        t_k_old = t_k;
        t_k = t_k_plus_1;
        
        x_k_old = x_k;
        x_k = x_k_plus_1;
        
        % --- 7. Update Total Work & Record History ---
        total_equivalent_work = total_equivalent_work + work_this_iter;
        if is_fine_level
            history.cumulative_work(k) = total_equivalent_work;
            history.iter(k) = k;
            history.obj_val(k) = obj_f_k(x_k) + obj_g_k(x_k);
            history.gt_error(k) = norm(x_k - x_star);
            
            % Check for convergence
            if k > 1
                change = norm(x_k - x_k_old) / norm(x_k);
                if change < params.fine_tolerance
                    fprintf('Convergence reached at fine-level iteration %d.\n', k);
                    history.iter = history.iter(1:k);
                    history.obj_val = history.obj_val(1:k);
                    history.gt_error = history.gt_error(1:k);
                    history.cumulative_work = history.cumulative_work(1:k);
                    break;
                end
            end
        else
            % Coarse-level convergence
            if k > 1 && (norm(x_k - x_k_old) / norm(x_k) < params.tolerance)
                break;
            end
        end
    end
    
    % Store total work for this level's call (for recursion)
    if ~is_fine_level
        history.cumulative_work(1) = total_equivalent_work;
    end
    
end