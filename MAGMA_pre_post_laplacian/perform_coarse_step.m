function [y_k_plus_1, s_k, total_equivalent_work, d_k] = perform_coarse_step(x_k, L_level, R_level, P_level, k, params, Lf_level, level_no, grad_F_mu_x_k, obj_f_x_k, v_H_parent)
global N_vec
%PERFORM_COARSE_STEP Implements a recursive V-cycle that carries the v_H term down.
%   The backtracking line search is not performed at any level of the recursion.

    % --- 1. Parameter and Operator Setup for the NEXT level down ---
    next_level = level_no + 1;
    L_fine   = L_level{level_no};
    L_coarse = L_level{next_level};
    R = R_level{level_no};
    P = P_level{level_no};
    lambda = params.lambda; mu = params.mu;
    NH = params.NH; c = params.c; tau = params.tau;

    % --- 2. Define the Coarse Problem for the NEXT level ---
    % Handle the initial call where no parent v_H exists.
    if isempty(v_H_parent)
        v_H_parent = zeros(size(x_k));
    end

    x_H0 = R * x_k;
    grad_F_mu_fine = grad_F_mu_x_k;
    
    % Calculate the initial gradient on the coarse grid
    grad_F_mu_coarse_initial = (1/N_vec(level_no)^2)*L_coarse' * (L_coarse * x_H0) + lambda * (x_H0 ./ sqrt(mu^2 + x_H0.^2));
    
    % Calculate the local v_H term for this level
    v_H_local = R * grad_F_mu_fine - grad_F_mu_coarse_initial;

    % **KEY CHANGE**: The total v_H is the local term plus the restricted parent term.
    v_H_total = v_H_local + R * v_H_parent;

    % Define the coarse objective and gradient using the total v_H term.
    obj_f_coarse = @(x_H) (1/N_vec(level_no)^2)*0.5 * norm(L_coarse * x_H)^2 + lambda * sum(sqrt(mu^2 + x_H.^2)) + dot(v_H_total, x_H);
    grad_f_coarse = @(x_H) (1/N_vec(level_no)^2)*L_coarse' * (L_coarse * x_H)+ lambda * (x_H ./ sqrt(mu^2 + x_H.^2)) + v_H_total;

    % --- Gradient Checking ---
    DEBUG_GRADIENT_CHECK = false; % Set to true to enable the check
    if DEBUG_GRADIENT_CHECK
        fprintf('--- Running Gradient Check for Level %d -> %d ---\n', level_no, next_level);
        check_point = x_H0;
        
        % Using the Optimization Toolbox's checkGradients correctly
        % This requires the helper function `gradient_checker_wrapper` below.
        [valid, err] = checkGradients(@(x) gradient_checker_wrapper(x, obj_f_coarse, grad_f_coarse), check_point);
        
        diff_norm = norm(err.Objective);
        fprintf('Norm of relative difference vector: %e\n', diff_norm);
        
        if ~valid
            warning('Gradient may be incorrect! checkGradients returned false.');
        else
            fprintf('Gradient check passed successfully.\n');
        end
        fprintf('-----------------------------------------------\n');
    end


    % --- 3. Recursive Solver Logic ---
    if next_level == length(L_level)
        % --- BASE CASE: We have reached the COARSEST level ---
        fprintf('    -> [Level %d] Solving at COARSEST level with MFISTA...\n', next_level);
        % fprintf('    -> [Level %d] Solving at COARSEST level with FminUnc...\n', next_level);
        prox_op_coarse = @(x, l) x;
        
        x_H_final = mfista_solver(obj_f_coarse, @(x) 0, grad_f_coarse, prox_op_coarse, x_H0, Lf_level{next_level}, lambda, NH,params.tolerance);
        % x_H_final = fminunc(obj_f_coarse,x_H0);
        coarsest_level_factor = get_work_factor(next_level);
        total_equivalent_work = NH * coarsest_level_factor;
    else
        % --- RECURSIVE STEP: Go one level deeper ---
        fprintf('    -> [Level %d] Recursing to solve problem at Level %d...\n', level_no, next_level);
        
        %pre-smoothing
        prox_op_coarse = @(x, l) x;
        x_H0_smooth = mfista_solver(obj_f_coarse, @(x) 0, grad_f_coarse, prox_op_coarse, x_H0, Lf_level{next_level}, lambda, NH,params.tolerance);

        % Also pass the full gradient of the new coarse objective.
        [x_H_final, s_k, work_from_below, d_k] = perform_coarse_step( ...
            x_H0_smooth, L_level, R_level, P_level, k, params, Lf_level, ...
            next_level, grad_f_coarse(x_H0_smooth), obj_f_coarse(x_H0_smooth), v_H_total ...
        );

        %post-smoothing
        prox_op_coarse = @(x, l) x;
        x_H_final = mfista_solver(obj_f_coarse, @(x) 0, grad_f_coarse, prox_op_coarse, x_H_final, Lf_level{next_level}, lambda, NH,params.tolerance);

        smooth_level_factor = get_work_factor(next_level);
        smooth_equivalent_work = NH * smooth_level_factor;

        total_equivalent_work = work_from_below + 2*smooth_equivalent_work;
        
    end

    %The correction is calculated and applied at each level on the way back up the recursion.

        d_k = P * (x_H_final - x_H0);
        
        % if level_no == 1
            % --- TOP LEVEL CALL: Perform the line search to find the optimal step size s_k ---
            fprintf('    -> [Level %d] Performing final backtracking line search...\n', level_no);
            s_k = 1.0; 
            F_mu_k = obj_f_x_k + lambda * sum(sqrt(mu^2 + x_k.^2)); 
            descent_term = dot(d_k, grad_F_mu_fine);
    
            if descent_term >= 0
                fprintf('    -> Warning: Coarse direction is not a descent direction. Taking zero step.\n');
                descent_term
                y_k_plus_1 = x_k; s_k = 0;
                return;
            else
                while true
                    
                    y_k_plus_1 = x_k + s_k * d_k;
                    F_mu_candidate = calculate_smoothed_objective(y_k_plus_1, L_fine, lambda, params);
                    if F_mu_candidate <= F_mu_k + c * s_k * descent_term
                        fprintf('    -> Line search success. Step size s_k = %.4f\n', s_k);
                        break; 
                    end
                    s_k = s_k * tau;
                    if s_k < 1e-4
                        fprintf('    -> Warning: Line search failed. Taking zero step.\n');
                        y_k_plus_1 = x_k; s_k = 0;
                        % y_k_plus_1 = max(0, y_k_plus_1); % Ensure non-negativity
                        return;
                    end
                end
            end
            % y_k_plus_1 = x_k + s_k * d_k; % Apply correction with optimal step size
            % Ensure non-negativity before returning
    
            % y_k_plus_1 = max(0, y_k_plus_1);
            
        % else
            % --- RECURSIVE SUB-CALL: Apply correction with a fixed step size of 1 ---
            % We don't do a line search here. We compute the corrected solution on this
            % grid and return it to the level above.
            % s_k = 1.0; % Use a fixed step size for intermediate levels
            % y_k_plus_1 = x_k + s_k*d_k; % Apply the full correction
        % end

end

    
%     if level_no == 1
%         % The correction is calculated and applied at each level on the way back up the recursion.
%         d_k = P * (x_H_final - x_H0);
% 
%         % --- TOP LEVEL CALL: Perform the line search to find the optimal step size s_k ---
%         s_k = 1.0; 
%         descent_term = dot(d_k, grad_F_mu_fine);
%         if descent_term >= 0
%             fprintf('    -> Warning: Coarse direction is not a descent direction. Taking zero step.\n');
%             y_k_plus_1 = x_k; s_k = 0;
%             return;
%         end
%         y_k_plus_1 = x_k + s_k * d_k; % Apply correction with step size 1
%         y_k_plus_1 = max(0, y_k_plus_1); % Ensure non-negativity
%     else
%         s_k = 1.0;
%         d_k = P * (x_H_final - x_H0);
%         y_k_plus_1 = x_k + s_k * d_k; % Apply correction with step size 1
%         y_k_plus_1 = max(0, y_k_plus_1); % Ensure non-negativity
%     end
% 
% end