% run_magma.m - Main Driver Script for MAGMA Tomography Reconstruction

clear; clc; close all;

echo off

%% 1. Problem Setup
sampleName = 'Phantom'; % choose from {'Phantom', 'Brain', ''Golosio', 'circle', 'checkboard', 'fakeRod'}; not case-sensitive

Nx = 127; Ny = Nx; WSz = [Ny, Nx]; %  XH: Ny = 100 -> 10 or 50 for debug.  currently assume object shape is square not general rectangle

WGT = createObject(sampleName, WSz); 

WGT = WGT/max(WGT(:)); assert(all(WGT(:)>=0), 'Groundtruth object should be nonnegative');% Object without noise. always assume WGT0 is normailzed with maximum 1.

NTheta = 20;  % 30/45/60/90 sample angle #. Use odd NOT even, for display purpose of sinagram of Phantom. As even NTheta will include theta of 90 degree where sinagram will be very bright as Phantom sample has verical bright line on left and right boundary.
NTau = ceil(sqrt(sum(WSz.^2))); NTau = NTau + rem(NTau-Ny,2); % number of discrete beam, enough to cover object diagonal, plus tolarence with maxDrift, also use + rem(NTau-Ny,2) to make NTau the same odd/even as Ny just for perfection, so that for theta=0, we have sum(WGT, 2)' and  S(1, (1:Ny)+(NTau-Ny)/2) are the same with a scale factor
SSz = [NTheta, NTau];

gaussianSTD = 0; % 0.01; % [0, 0.01, 0.05] % noise

L = XTM_Tensor_XH(WSz, NTheta, NTau, 0, WGT);
LNormalizer = full(max(sum(L,2)));
L = L/LNormalizer;

S = reshape(L*WGT(:), NTheta, NTau);
SMax = max(S(:));        
rng('default'); % same random number for initial test
S = imnoise(S/SMax, 'gaussian', 0, gaussianSTD^2)*SMax; % add noise to [0, 1 ] grayscale image, imnoise add gaussian noise but also seems to still threshold the noisy image to [0, 1]                 

figNo = 1;
figure(figNo+1); imagesc(WGT); axis image; title('Object Ground Truth');
figure(figNo+2); multAxes(@imagesc, S); title('Sinogram');linkAxesXYZLimColorView(); 

s = reshape(S,NTau*NTheta,1);

x_star = reshape(WGT, Nx*Nx,1); %ground truth
%% 2. Algorithm Parameters
% These are based on the paper's experiments (Section 4.3)
params.lambda = 1e-6;   % Regularization parameter
params.kappa = 0.1;     % Coarse direction selection threshold (3.4)
params.theta = 1;     % Proximity threshold (e.g., 0.1 * ||x||) (3.4)
params.Kd = 30;         % Max consecutive gradient steps (3.4)
params.NH = 10;         % Number of iterations for the coarse-level solver
params.T = 200;         % Max iterations for the main MAGMA loop
params.mu = 1e-4;       % NEW: Smoothing parameter for the gradient calculation
params.c = 0.1;         % Armijo line search parameter, a typical value
params.tau = 0.5;           % Line search shrinkage factor (e.g., 0.95)
params.tolerance = 1e-3; % Convergence tolerance for ||D(x)||
params.fine_tolerance = 1e-4;
% Each recursive level will perform 50% more iterations than the level above it.
% Level 1 runs for T=50.
% Level 2 will run for T=floor(50 * 1.5) = 75.
% Level 3 will run for T=floor(75 * 1.5) = 112, and so on.
params.magma_factor = 1; 
params.fista_factor = 1; 

%% 3. Multilevel Operator Setup
fprintf('Setting up multilevel operators...\n');
num_levels = 2; % Example: 128 -> 64 -> 32 -> 16 -> 8
[L_level, R_level, P_level, N_vec, Lf_level] = setup_multilevel_operators_new(L, Nx, num_levels);
% L_level{1} is your fine matrix L
% R_level{1} is the restriction op from level 1 to 2
% P_level{1} is the prolongation op from level 2 to

%% 4. Initialization
fprintf('Initializing variables...\n');
x0 = zeros(Nx^2, 1); % Initial guess is a black image

%% 5. Run the MAGMA Solver
fprintf('Starting MAGMA solver...\n');

% In main script...

% --- Define Fine-Level Functions ---
obj_f_fine  = @(x) 0.5 * norm(L_level{1} * x - s)^2;
obj_g_fine  = @(x) params.lambda * norm(x, 1);
grad_f_fine = @(x) L_level{1}' * (L_level{1} * x - s);
prox_op_fine = @prox_g;

% --- Add full operator lists to params for easy access in recursion ---
params.L_level_full = L_level;
params.R_level_full = R_level;
params.P_level_full = P_level;
params.Lf_level_full = Lf_level;

% Run the MAGMA Solver
level_no = 1;% strating with level no 1
fprintf('Starting MAGMA solver...\n');
tic;

[x_recon_magma, history_magma] = magma_solver(obj_f_fine, obj_g_fine, grad_f_fine, prox_op_fine, ...
                                             L_level{1}, Lf_level{1}, R_level{1}, ...
                                             s, x0, params, level_no, x_star);
Magma_time = toc;

fprintf('MAGMA reconstruction complete in %.2f s.\n', Magma_time);

%%
% Run FISTA Solver
max_iters = params.T;
obj_f_fine  = @(x) 0.5 * norm(L_level{1} * x - s)^2;
obj_g_fine  = @(x) params.lambda * norm(x, 1);
grad_f_fine = @(x) L_level{1}' * (L_level{1} * x - s);
prox_op_fine = @prox_g;
fprintf('Starting FISTA solver...\n');

[x_recon_fista, history_fista] = fista_solver_instrumented(obj_f_fine, obj_g_fine, grad_f_fine, @prox_g, x0,Lf_level{1}, ...
    params, max_iters, x_star);

Fista_time = history_fista.elapsed_time;
fprintf('FISTA reconstruction complete in %.2f s.\n', Fista_time);
%%
% --- Run AGM Solver ---
fprintf('Starting AGM solver...\n');
tic;

[x_recon_agm, history_agm] = agm_solver(obj_f_fine, obj_g_fine, grad_f_fine, @prox_g, ...
                                        Lf_level{1}, s, x0, params, x_star);
AGM_time = toc;
fprintf('AGM reconstruction complete in %.2f s.\n', AGM_time);


%% 6. Display Results and Comparison Metrics

fprintf('\n--- 6. Calculating final metrics and displaying results ---\n');

% --- Reshape final vectors into images ---
W_magma_recon = reshape(x_recon_magma, Ny, Nx);
W_fista_recon = reshape(x_recon_fista, Ny, Nx);
W_agm_recon = reshape(x_recon_agm, Ny, Nx);

% --- Calculate PSNR values ---

psnr_magma = psnr(W_magma_recon, WGT);
psnr_fista = psnr(W_fista_recon, WGT);
psnr_agm = psnr(W_agm_recon, WGT);

% --- Assume computation times are in these variables from your main script ---
% Magma_time should be from the external tic/toc in your main script.

% --- Figure 1: Visual Comparison of Reconstructed Images ---
figure('Name', 'Image Reconstruction Comparison', 'Position', [100, 100, 1800, 500]);
subplot(1, 4, 1); imagesc(WGT); axis image; title('Ground Truth'); xlabel(sprintf('(%d x %d)', Nx, Ny));
subplot(1, 4, 2); imagesc(W_fista_recon); axis image; title(sprintf('FISTA\nPSNR: %.2f dB, Time: %.2f s', psnr_fista, Fista_time));
subplot(1, 4, 3); imagesc(W_agm_recon); axis image; title(sprintf('AGM\nPSNR: %.2f dB, Time: %.2f s', psnr_agm, AGM_time));
subplot(1, 4, 4); imagesc(W_magma_recon); axis image; title(sprintf('MAGMA\nPSNR: %.2f dB, Time: %.2f s', psnr_magma, Magma_time));

figure('Name', 'Algorithm Convergence Comparison', 'Position', [200, 200, 1400, 600]);
subplot(1, 2, 1);
semilogy(history_fista.iter, history_fista.obj_val, 'b-o', 'DisplayName', 'FISTA'); hold on;
semilogy(history_agm.iter, history_agm.obj_val, 'g-d', 'DisplayName', 'AGM (User)');
semilogy(history_magma.cumulative_work, history_magma.obj_val, 'r-s', 'LineWidth', 2, 'DisplayName', 'MAGMA');
hold off; grid on; title('Objective Function Convergence'); xlabel('Equivalent Fine-Level Work'); ylabel('Objective Value (log scale)'); legend;

subplot(1, 2, 2);
semilogy(history_fista.iter, history_fista.gt_error, 'b-o', 'DisplayName', 'FISTA'); hold on;
semilogy(history_agm.iter, history_agm.gt_error, 'g-d', 'DisplayName', 'AGM (User)');
semilogy(history_magma.cumulative_work, history_magma.gt_error, 'r-s', 'LineWidth', 2, 'DisplayName', 'MAGMA');
hold off; grid on; title('Ground Truth Error vs. Iteration Count'); xlabel('Equivalent Fine-Level Work'); ylabel('Ground Truth Error (log scale)'); legend;

%% 7. Update Step Visualization (5-Way Comparison with AGM)
fprintf('\n--- 7. Visualizing update steps ---\n');

% --- Select iterations to plot from the general iteration list ---

iters_to_plot = 1:min([length(history_magma.iter), length(history_fista.iter), length(history_agm.iter)]);
num_available = length(iters_to_plot);

if num_available == 0
    fprintf('No iterations available to plot.\n');
else
    % Select up to 5 iterations if available
    if num_available >= 5
        num_plots = 5;
        indices = round(linspace(1, num_available, num_plots));
    else
        num_plots = num_available;
        indices = 1:num_available;
    end
    
    iters_to_plot_final = iters_to_plot(indices);

    % --- Figure 3: Image Comparison (5-Way) ---
    figure('Name', 'Update Step Image Comparison (5-Way)', 'Position', [50, 50, 1500, 250 * num_plots]);
    sgtitle('Side-by-Side Image Comparison of Algorithm State Vectors', 'FontSize', 16, 'FontWeight', 'bold');
    
    for i = 1:length(iters_to_plot_final)
        iter_num = iters_to_plot_final(i);
        
        % Check if this iteration was a coarse step to get the correct MAGMA data index
        cell_idx = find(history_magma.coarse_step_iters == iter_num, 1);
        if isempty(cell_idx)
            fprintf('Info: Skipping plot for iteration %d because it was not a MAGMA coarse step.\n', iter_num);
            continue; % Go to the next iteration
        end

        % Retrieve all five vectors and reshape
        img_fista   = reshape(history_fista.y_k_hist(:, iter_num), Ny, Nx);
        img_agm     = reshape(history_agm.y_k_hist(:, iter_num), Ny, Nx);
        img_magma_y = reshape(history_magma.coarse_step_updates{cell_idx}, Ny, Nx);
        img_magma_x = reshape(history_magma.x_k_coarse_hist{cell_idx}, Ny, Nx);
        img_magma_z = reshape(history_magma.z_k_coarse_hist{cell_idx}, Ny, Nx);

        % Correct 5-way subplot syntax
        subplot(num_plots, 5, 5*i - 4); imagesc(img_fista); axis image; title(sprintf('FISTA y_k (Iter: %d)', iter_num));
        subplot(num_plots, 5, 5*i - 3); imagesc(img_agm); axis image; title('AGM y_k');
        subplot(num_plots, 5, 5*i - 2); imagesc(img_magma_y); axis image; title('MAGMA y_{k+1}');
        subplot(num_plots, 5, 5*i - 1); imagesc(img_magma_x); axis image; title('MAGMA x_k');
        subplot(num_plots, 5, 5*i - 0); imagesc(img_magma_z); axis image; title('MAGMA z_{k+1}');
    end

    % --- Figure 4: Pixel-wise Plot (5-Way) ---
    figure('Name', 'Update Step Pixel Plot (5-Way)', 'Position', [100, 100, 2000, 250 * num_plots]);
    sgtitle('Side-by-Side Pixel-wise Comparison of Algorithm State Vectors', 'FontSize', 16, 'FontWeight', 'bold');
    
    for i = 1:length(iters_to_plot_final)
        iter_num = iters_to_plot_final(i);

        % Repeat the check for the second plot
        cell_idx = find(history_magma.coarse_step_iters == iter_num, 1);
        if isempty(cell_idx)
            continue;  
        end

        % Retrieve all five vectors
        y_fista = history_fista.y_k_hist(:, iter_num);
        y_agm   = history_agm.y_k_hist(:, iter_num);
        y_magma = history_magma.coarse_step_updates{cell_idx};
        x_magma = history_magma.x_k_coarse_hist{cell_idx};
        z_magma = history_magma.z_k_coarse_hist{cell_idx};

        % Correct 5-way subplot syntax
        subplot(num_plots, 5, 5*i - 4); plot(y_fista, 'b-'); grid on; xlim([1, length(y_fista)]); title(sprintf('FISTA y_k (Iter: %d)', iter_num)); ylabel('Value');
        subplot(num_plots, 5, 5*i - 3); plot(y_agm, 'g-'); grid on; xlim([1, length(y_agm)]); title('AGM y_k');
        subplot(num_plots, 5, 5*i - 2); plot(y_magma, 'r-'); grid on; xlim([1, length(y_magma)]); title('MAGMA y_{k+1} [Coarse]');
        subplot(num_plots, 5, 5*i - 1); plot(x_magma, 'c-'); grid on; xlim([1, length(x_magma)]); title('MAGMA x_k');
        subplot(num_plots, 5, 5*i - 0); plot(z_magma, 'm-'); grid on; xlim([1, length(z_magma)]); title('MAGMA z_{k+1}');
        
        if i == length(iters_to_plot_final)
             subplot(num_plots, 5, 5*i - 4); xlabel('Pixel Index');
             subplot(num_plots, 5, 5*i - 3); xlabel('Pixel Index');
             subplot(num_plots, 5, 5*i - 2); xlabel('Pixel Index');
             subplot(num_plots, 5, 5*i - 1); xlabel('Pixel Index');
             subplot(num_plots, 5, 5*i - 0); xlabel('Pixel Index');
        end
    end
end
%%

% %% 7. Update Step Visualization (4-Way Side-by-Side Layout)
% fprintf('\n--- 7. Visualizing update steps (4-Way Side-by-Side Layout) ---\n');
% 
% % --- Select iterations to plot ---
% iters_to_plot = 1:min(length(history_magma.iter),length(history_fista.iter)) -1;
% num_available = length(iters_to_plot);
% 
% if num_available == 0
%     fprintf('No coarse steps were taken by MAGMA, skipping update visualization.\n');
% else
%     % Select up to 5 iterations if available
%     if num_available >= 5
%         num_plots = 5;
%         indices = round(linspace(1, num_available, num_plots));
%     else
%         num_plots = num_available;
%         indices = 1:num_available;
%     end
% 
%     iters_to_plot = iters_to_plot(indices);
%     selected_cell_indices = indices; 
% 
%     % --- Figure 3: Image Comparison of Update Vectors (Nx4 Grid) ---
%     figure('Name', 'Update Step Image Comparison (4-Way)', 'Position', [50, 50, 1200, 250 * num_plots]);
%     sgtitle('Side-by-Side Image Comparison of Algorithm State Vectors', 'FontSize', 16, 'FontWeight', 'bold');
% 
%     for i = 1:num_plots
%         iter_num = iters_to_plot(i);
%         cell_idx = selected_cell_indices(i);
% 
%         % Retrieve all four vectors and reshape to images
%         img_fista_y = reshape(history_fista.y_k_hist(:, iter_num), Ny, Nx);
%         img_magma_y = reshape(history_magma.coarse_step_updates{cell_idx}, Ny, Nx);
%         img_magma_x = reshape(history_magma.x_k_coarse_hist{cell_idx}, Ny, Nx);
%         img_magma_z = reshape(history_magma.z_k_coarse_hist{cell_idx}, Ny, Nx);
% 
%         subplot(num_plots, 4, 4*i - 3); imagesc(img_fista_y); axis image; colormap gray; title(sprintf('FISTA y_k (Iter: %d)', iter_num));
%         subplot(num_plots, 4, 4*i - 2); imagesc(img_magma_y); axis image; colormap gray; title(sprintf('MAGMA y_{k+1} [Coarse]'));
%         subplot(num_plots, 4, 4*i - 1); imagesc(img_magma_x); axis image; colormap gray; title(sprintf('MAGMA x_k'));
%         subplot(num_plots, 4, 4*i - 0); imagesc(img_magma_z); axis image; colormap gray; title(sprintf('MAGMA z_{k+1}'));
%     end
% 
%     % --- Figure 4: Pixel-wise Plot of Update Vectors (Nx4 Grid) ---
%     figure('Name', 'Update Step Pixel Plot (4-Way)', 'Position', [100, 100, 1600, 250 * num_plots]);
%     sgtitle('Side-by-Side Pixel-wise Comparison of Algorithm State Vectors', 'FontSize', 16, 'FontWeight', 'bold');
% 
%     for i = 1:num_plots
%         iter_num = iters_to_plot(i);
%         cell_idx = selected_cell_indices(i);
% 
%         % Retrieve all four vectors
%         y_fista = history_fista.y_k_hist(:, iter_num);
%         y_magma = history_magma.coarse_step_updates{cell_idx};
%         x_magma = history_magma.x_k_coarse_hist{cell_idx};
%         z_magma = history_magma.z_k_coarse_hist{cell_idx};
% 
%         % Plot FISTA y_k
%         subplot(num_plots, 4, 4*i - 3);
%         plot(y_fista, 'b-'); grid on; xlim([1, length(y_fista)]);
%         title(sprintf('FISTA y_k (Iter: %d)', iter_num)); ylabel('Value');
% 
%         % Plot MAGMA y_k+1
%         subplot(num_plots, 4, 4*i - 2);
%         plot(y_magma, 'r-'); grid on; xlim([1, length(y_magma)]);
%         title('MAGMA y_{k+1} [Coarse]');
% 
%         % Plot MAGMA x_k
%         subplot(num_plots, 4, 4*i - 1);
%         plot(x_magma, 'g-'); grid on; xlim([1, length(x_magma)]);
%         title('MAGMA x_k');
% 
%         % Plot MAGMA z_k+1
%         subplot(num_plots, 4, 4*i - 0);
%         plot(z_magma, 'm-'); grid on; xlim([1, length(z_magma)]);
%         title('MAGMA z_{k+1}');
% 
%         % Add xlabel to bottom row only
%         if i == num_plots
%             subplot(num_plots, 4, 4*i - 3); xlabel('Pixel Index');
%             subplot(num_plots, 4, 4*i - 2); xlabel('Pixel Index');
%             subplot(num_plots, 4, 4*i - 1); xlabel('Pixel Index');
%             subplot(num_plots, 4, 4*i - 0); xlabel('Pixel Index');
%         end
%     end
% end