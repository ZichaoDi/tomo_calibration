% run_mr_fista_comparison.m - Main Driver Script for MR-FISTA Tomography
clear; clc; close all;
echo off

%% 1. Problem Setup
sampleName = 'Phantom'; % choose from {'Phantom', 'Brain', ''Golosio', 'circle', 'checkboard', 'fakeRod'}; not case-sensitive

Nx = 256; Ny = Nx; WSz = [Ny, Nx]; %  XH: Ny = 100 -> 10 or 50 for debug.  currently assume object shape is square not general rectangle

WGT = createObject(sampleName, WSz); 

WGT = WGT/max(WGT(:)); assert(all(WGT(:)>=0), 'Groundtruth object should be nonnegative');% Object without noise. always assume WGT0 is normailzed with maximum 1.

NTheta = 50;  % 30/45/60/90 sample angle #. Use odd NOT even, for display purpose of sinagram of Phantom. As even NTheta will include theta of 90 degree where sinagram will be very bright as Phantom sample has verical bright line on left and right boundary.
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
params.lambda = 1e-4;   % Regularization parameter (higher for noisy data)
params.kappa = 0.1;     % Coarse direction selection threshold (from paper)
params.theta = 1e-3;    % Proximity threshold (relative, from SolveMRFISTA)
params.Kd = 30;         % Max consecutive gradient steps (from paper)
params.NH = 10;         % Number of iterations for coarse-level solver
params.T = 200;         % Max iterations for the main loop
params.mu= 1e-4;        % Smoothing parameter for gradient
params.c = 1e-4;        % Armijo line search parameter (a typical value)
params.tau = 0.5;       % Line search shrinkage factor
params.tolerance = 1e-6; % Convergence tolerance for coarse levels
params.fine_tolerance = 1e-6; % Convergence for fine level

%% 3. Multilevel Operator Setup
fprintf('Setting up multilevel operators...\n');
num_levels = 3; % 128 -> 64 -> 32
[L_level, R_level, P_level, N_vec, Lf_level] = setup_multilevel_operators_new(L, Nx, num_levels);
params.N_vec = N_vec; % Add N_vec to params for work calculation

%% 4. Initialization
fprintf('Initializing variables...\n');
x0 = zeros(Nx^2, 1); % Initial guess is a black image

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
params.s = s; % Add sinogram to params

%% 5. Run the MR-FISTA Solver
level_no = 1; % Start at the finest level (level 1)
v0 = zeros(Nx^2, 1); % Initial coherency term is zero
fprintf('Starting MR-FISTA solver...\n');
tic;
[x_recon_mrfista, history_mrfista] = mrfista_solver(obj_f_fine, obj_g_fine, grad_f_fine, prox_op_fine, ...
                                                    x0, params, level_no, x_star, v0);
MRFISTA_time = toc;
fprintf('MR-FISTA reconstruction complete in %.2f s.\n', MRFISTA_time);

%% 6. Run FISTA Solver for Comparison
max_iters = params.T * params.fista_factor; % Equivalent work
fprintf('Starting FISTA solver...\n');
tic;
[x_recon_fista, history_fista] = fista_solver_instrumented(obj_f_fine, obj_g_fine, grad_f_fine, @prox_g, x0, Lf_level{1}, ...
    params, max_iters, x_star);
Fista_time = toc;
fprintf('FISTA reconstruction complete in %.2f s.\n', Fista_time);

%% 7. Run AGM Solver for Comparison
fprintf('Starting AGM solver...\n');
max_iters_agm = params.T * params.fista_factor;
tic;
[x_recon_agm, history_agm] = agm_solver(obj_f_fine, obj_g_fine, grad_f_fine, @prox_g, ...
                                        Lf_level{1}, s, x0, params, x_star, max_iters_agm);
AGM_time = toc;
fprintf('AGM reconstruction complete in %.2f s.\n', AGM_time);

%% 8. Display Results and Comparison Metrics
fprintf('\n--- 8. Calculating final metrics and displaying results ---\n');

% --- Reshape final vectors into images ---
W_mrfista_recon = reshape(x_recon_mrfista, Ny, Nx);
W_fista_recon = reshape(x_recon_fista, Ny, Nx);
W_agm_recon = reshape(x_recon_agm, Ny, Nx);

% --- Calculate PSNR values ---
psnr_mrfista = psnr(W_mrfista_recon, WGT);
psnr_fista = psnr(W_fista_recon, WGT);
psnr_agm = psnr(W_agm_recon, WGT);

% --- Figure 1: Visual Comparison of Reconstructed Images ---
figure('Name', 'Image Reconstruction Comparison', 'Position', [100, 100, 1800, 500]);
subplot(1, 4, 1); imagesc(WGT); axis image; colormap gray; title('Ground Truth'); xlabel(sprintf('(%d x %d)', Nx, Ny));
subplot(1, 4, 2); imagesc(W_fista_recon); axis image; colormap gray; title(sprintf('FISTA\nPSNR: %.2f dB, Time: %.2f s', psnr_fista, Fista_time));
subplot(1, 4, 3); imagesc(W_agm_recon); axis image; colormap gray; title(sprintf('AGM\nPSNR: %.2f dB, Time: %.2f s', psnr_agm, AGM_time));
subplot(1, 4, 4); imagesc(W_mrfista_recon); axis image; colormap gray; title(sprintf('MR-FISTA\nPSNR: %.2f dB, Time: %.2f s', psnr_mrfista, MRFISTA_time));

% --- Figure 2: Convergence Plots ---
figure('Name', 'Algorithm Convergence Comparison', 'Position', [200, 200, 1400, 600]);
subplot(1, 2, 1);
semilogy(history_fista.iter, history_fista.obj_val, 'b-o', 'DisplayName', 'FISTA', 'MarkerSize', 4, 'LineWidth', 1); hold on;
semilogy(history_agm.iter, history_agm.obj_val, 'g-d', 'DisplayName', 'AGM', 'MarkerSize', 4, 'LineWidth', 1);
semilogy(history_mrfista.cumulative_work, history_mrfista.obj_val, 'r-s', 'LineWidth', 2, 'DisplayName', 'MR-FISTA');
hold off; grid on; title('Objective Function Convergence'); xlabel('Equivalent Fine-Level Work'); ylabel('Objective Value (log scale)'); legend;
ylim([min(history_mrfista.obj_val)*0.9, history_fista.obj_val(1)*1.1]);

subplot(1, 2, 2);
semilogy(history_fista.iter, history_fista.gt_error, 'b-o', 'DisplayName', 'FISTA', 'MarkerSize', 4, 'LineWidth', 1); hold on;
semilogy(history_agm.iter, history_agm.gt_error, 'g-d', 'DisplayName', 'AGM', 'MarkerSize', 4, 'LineWidth', 1);
semilogy(history_mrfista.cumulative_work, history_mrfista.gt_error, 'r-s', 'LineWidth', 2, 'DisplayName', 'MR-FISTA');
hold off; grid on; title('Ground Truth Error vs. Work'); xlabel('Equivalent Fine-Level Work'); ylabel('||x_k - x*||_2 (log scale)'); legend;
ylim([min(history_mrfista.gt_error)*0.9, history_fista.gt_error(1)*1.1]);