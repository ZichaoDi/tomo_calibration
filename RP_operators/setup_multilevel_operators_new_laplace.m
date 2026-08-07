function [L_level, R_level, P_level, N_vec, Lf_level] = setup_multilevel_operators_new(L_fine, Nx, num_levels)
global N_vec
% Purpose:
%   Generates system matrices (L), restriction (R), and prolongation (P)
%   operators.
%   --- NEW ---
%   Also pre-calculates the Lipschitz constant (Lf) and the L'*L matrix
%   for each level to improve solver performance.
%
% Inputs:
%   L_fine: The system matrix for the finest grid level.
%   Nx: The dimension of the finest grid image (assuming Nx = Ny).
%   num_levels: The total number of levels in the hierarchy.
%
% Outputs:
%   L_level: Cell array, L_level{i} is the system matrix for level i.
%   R_level: Cell array, R_level{i} restricts from level i to i+1.
%   P_level: Cell array, P_level{i} prolongates from level i+1 to i.
%   N_vec:   Vector of grid sizes for each level.
%   Lf_level:  --- NEW --- Cell array for the Lipschitz constant for each level.
%   LtL_level: --- NEW --- Cell array for the L'*L matrix for each level.

    fprintf('--- Setting up %d levels of multigrid operators and constants ---\n', num_levels);

    global N;
    N = zeros(1, num_levels);
    N(1) = Nx;
    for i = 2:num_levels
        N(i) = floor((N(i-1) + 1) / 2);
    end
    N_vec = N;
    fprintf('Grid hierarchy (N): %s\n', num2str(N));

    % Initialize output
    L_level = cell(num_levels, 1);
    R_level = cell(num_levels - 1, 1);
    P_level = cell(num_levels - 1, 1);
    Lf_level = cell(num_levels, 1);
    % LtL_level = cell(num_levels, 1);

    % -- Make finest level sparse --
    L_level{1} = sparse(L_fine);
    % LtL_level{1} = sparse(L_level{1}' * L_level{1});
    Lf_level{1} = estimate_lipschitz_power(L_level{1});
    fprintf('  -> Lf for level 1 set to: %.2f\n', Lf_level{1});

    for level = 1:(num_levels - 1)
        fprintf('Generating operators and constants for level %d -> %d\n', level, level + 1);

        m_fine = N(level);
        IhH_1D = sparse(restrict_operator(m_fine, 'standard_FW'));

        R = sparse(kron(IhH_1D, IhH_1D));
        P = sparse(R');

        R_level{level} = (1/4)*R;
        P_level{level} = (1)*P;
        fprintf('  -> Stored R and P for level %d.\n', level);

        L_fine_current = L_level{level};
        L_coarse_next = sparse(R*L_fine_current * R');
        L_level{level + 1} = L_coarse_next;
        fprintf('  -> Generated L for level %d.\n', level + 1);

        % LtL_level{level + 1} = sparse(L_coarse_next' * L_coarse_next);
        Lf_level{level + 1} = estimate_lipschitz_power(L_coarse_next);
        fprintf('  -> Lf for level %d set to: %.2f\n', level + 1, Lf_level{level + 1});
    end

    fprintf('--- Operator and constant setup complete. ---\n');
    clear global N;
end