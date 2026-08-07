function [x, history, equivalent_work_done] = magma_solver(obj_f, obj_g, grad_f, prox_op, L_current, Lf_current, R_current, s, x0, params, level_no, x_star)
%MAGMA_SOLVER General MAGMA solver engine that tracks equivalent work.
    
    % --- Unpack Parameters ---
    lambda = params.lambda; kappa = params.kappa; theta = params.theta;
    Kd = params.Kd; T = params.T; NH = params.NH;
    tolerance = params.tolerance;
    is_fine_level = (level_no == 1);
    mu = params.mu;
    
    % Get the work factor for the current level
    my_level_factor = get_work_factor(level_no);
    
    % --- Initialization ---
    x_k = x0;
    y_k = x0; z_k = x0;
    eta_k = Lf_current;
    alpha_k = 0; 
    q = 0; x_tilde = x0;
    
    % --- History & Work Allocation ---
    history = [];
    total_equivalent_work = 0;
    if is_fine_level
        history.iter = zeros(T, 1);
        history.obj_val = zeros(T, 1);
        history.gt_error = zeros(T, 1);
        history.x_change = zeros(T, 1);
        history.cumulative_work = zeros(T, 1);

        history.coarse_step_updates = {}; % For y_k+1 from coarse step
        history.coarse_step_iters = [];   % Iteration numbers for coarse steps
        history.x_k_coarse_hist = {}; %  (For x_k)
        history.z_k_coarse_hist = {}; %  (For z_k+1)

        fprintf('Starting MAGMA on Fine Level 1 (T=%d)...\n', T);
    end

    % --- Main MAGMA Loop ---
    for k = 1:T
        x_k_old = x_k;
        
        if alpha_k == 0, t_k = 1; else, t_k = 1 / (alpha_k * eta_k); end
        t_k = min(max(t_k, 0), 1);
        fprintf('    -> t_k = %d.\n',  t_k);
        x_k = (1 - t_k) * y_k + t_k * z_k;
        
        grad_at_xk = grad_f(x_k);
        
        % for 3 or more level; the objective funciotn is already smooth so no need smooth part
        grad_g_mu = (lambda * x_k) ./ sqrt(mu^2 + x_k.^2);
        

        grad_F_mu =  grad_at_xk + grad_g_mu;

        obj_f_x_k = obj_f(x_k);
        
        % Decision Block
        cond1 = norm(R_current * grad_F_mu) > kappa * norm(grad_F_mu);
        cond2 = (norm(x_k - x_tilde) > theta * norm(x_tilde)) || (q < Kd);
        
        if cond1 && cond2 && level_no < length(params.L_level_full)
            % --- Coarse Correction Step ---
            [y_k_plus_1, s_k, coarse_work,d_k] = perform_coarse_step(x_k, ...
                params.L_level_full, s, params.R_level_full, params.P_level_full, k, params, ...
                params.Lf_level_full, level_no, grad_F_mu, obj_f_x_k);
            % Store the relevant vectors from this step
            if is_fine_level
                history.coarse_step_updates{end+1} = d_k;
                history.x_k_coarse_hist{end+1} = x_k; % <<< ADD THIS LINE
                history.coarse_step_iters(end+1) = k;
            end

            q = 0; x_tilde = x_k;
            % 
            eta_k_plus_1 = Lf_current;
            alpha_k_plus_1 = (k + 2) / (2 * Lf_current);

            % if s_k == 0
            %     eta_k_plus_1 = eta_k;
            %     alpha_k_plus_1 = alpha_k;
            % else
            %     if alpha_k == 0 % Should not happen due to "k > 1" but good for safety
            %         eta_k_plus_1 = params.Lf_level_full{level_no+1} / (params.c * s_k * kappa^2);
            %     else
            %         eta_k_plus_1 = max(1 / (4 * alpha_k^2 * eta_k), params.Lf_level_full{level_no+1} / (params.c * s_k * kappa^2));
            %     end
            %     alpha_k_plus_1 = (1 + sqrt(1 + 4 * eta_k_plus_1 * alpha_k^2 * eta_k)) / (2 * eta_k_plus_1);
            % end

            
            
            % Post-refinement FISTA runs for NH iterations at the CURRENT level

            % 
            % if is_fine_level
            %     NH =5;
            % end
            % 
            % y_k_plus_1 = fista_solver(grad_f, prox_op, y_k_plus_1, Lf_current, lambda, 2);
            % refinement_work = my_level_factor * NH;
            

            % Total work for this iteration
            work_this_iter = coarse_work + my_level_factor* 1; %+ refinement_work + my_level_factor * 1;
        else
            % --- Gradient Step (Proximal or Standard) ---
            fprintf('    -> [Level %d] Perfroming prox grad step...\n', level_no);
            grad_at_xk = grad_f(x_k);
            y_k_plus_1 = prox_op(x_k - (1/Lf_current) * grad_at_xk, lambda/Lf_current);
            y_k_plus_1 = max(0, y_k_plus_1);
            q = q + 1;
            
            eta_k_plus_1 = Lf_current;
            alpha_k_plus_1 = (k + 2) / (2 * Lf_current);
            
            % Work for a simple gradient step is 1 iteration at the current level
            work_this_iter = my_level_factor * 1;
        end
        
        % --- Mirror Descent & State Updates ---
        % grad_for_mirror = grad_f(x_k);
        z_k_plus_1 = perform_mirror_step(z_k, x_k, alpha_k_plus_1,  grad_at_xk, prox_op, lambda, is_fine_level);

        % <<<THIS BLOCK TO STORE z_k+1 IF A COARSE STEP WAS TAKEN >>>
        if is_fine_level
            history.z_k_coarse_hist{end+1} = z_k_plus_1;
        end
        
        y_k = y_k_plus_1; z_k = z_k_plus_1;
        
        eta_k = eta_k_plus_1; alpha_k = alpha_k_plus_1;

        % fprintf('eta_k = %d. \n', eta_k);
        % fprintf('alpha_k = %d. \n', alpha_k);
        
        % --- Update Total Work & Record History ---
        total_equivalent_work = total_equivalent_work + work_this_iter;

        if is_fine_level
            history.cumulative_work(k) = total_equivalent_work;
            history.iter(k) = k;
            history.obj_val(k) = obj_f_x_k + obj_g(x_k);
            history.gt_error(k) = norm(x_k - x_star);
            history.x_change(k) = norm1(x_k - x_k_old) / numel(x_k);
            
            if history.x_change(k) < params.fine_tolerance && k > 1
                fprintf('Convergence reached at fine-level iteration %d.\n', k);
                history.iter = history.iter(1:k);
                history.obj_val = history.obj_val(1:k);
                history.gt_error = history.gt_error(1:k);
                history.x_change = history.x_change(1:k);
                history.cumulative_work = history.cumulative_work(1:k);
                break;
            end
        else
            if k > 1 && norm(x_k - x_k_old) / numel(x_k) < tolerance
                fprintf('Convergence reached at coarse-level iteration %d.\n', k);
                break;
            end
        end
    end
    
    x = y_k;
    equivalent_work_done = total_equivalent_work;
end