function IhH=restrict_operator(m,frequency_type)
global N
j=find(N==m);

% Determine the size of the fine grid (nf) and coarse grid (nc)
nf = N(j);    % Fine grid size
nc = N(j+1);  % Coarse grid size

% Initialize the restriction operator matrix
I_h_2h = zeros(nc, nf);

% These 'c' variables seem to be part of a custom weighting or
% a frequency-dependent aspect. Assuming c=1 for now based on your code.
% If 'abs_FW' changes 'c', it needs to be applied consistently.
if(j==1 && strcmp(frequency_type,'abs_FW'))
    c=1; % Or whatever custom value applies for this case
else
    c=1;
end

% The logic needs to be based on the relationship between nf and nc
% For a standard full-weighting restriction:
% Coarse point i maps to fine points 2*i-1, 2*i, 2*i+1
% with weights typically 1/4, 1/2, 1/4 for interior points.
% Boundary points have different stencils.

% Your original code's if/else on mod(N(1),2) is problematic.
% The type of restriction operator should be chosen based on the
% grid point relationship, not just the very first grid size.

% Let's assume you intend a standard full-weighting restriction (1/4, 1/2, 1/4)
% for general cases where the fine grid is approximately twice the coarse grid.
% The indices used for `I_h_2h(i, ...)` should reflect the mapping.

% Common scenario for full-weighting (2h from h):
% Coarse point i corresponds to fine point 2i.
% Value at coarse point i comes from:
% 1/4 * value(fine_2i-1) + 1/2 * value(fine_2i) + 1/4 * value(fine_2i+1)

% Let's correct the loop for the full-weighting case.
% This implies nf should be approximately 2*nc - 1 (for an odd number of fine points)
% or 2*nc + 1 (for an odd number of coarse points).
% Or simply, if you are always halving and adding 1 (e.g., 127 -> 64),
% then nf = 2*nc - 1.

% Simplified and corrected logic for full-weighting (assuming 1/4, 1/2, 1/4 stencil)
% You had 1/2, 1, 1/2, which is another common stencil but needs careful application.
% Let's use your (1/2, 1, 1/2) stencil as it was in your code, but fix indexing.
% This stencil is often applied as:
% u_coarse(i) = 0.5*u_fine(2i-1) + 1.0*u_fine(2i) + 0.5*u_fine(2i+1)
% This requires nf to be 2*nc + 1, or special handling for boundaries if nf = 2*nc.

% Let's assume your intention was for the mapping to be:
% Coarse point i comes from fine points around 2*i.

% The problem statement implies N(1)=127, N(2)=64, N(3)=32.
% This means N(j) is the fine grid size, N(j+1) is the coarse grid size.

% Case 1: From 127 to 64 (nf=127, nc=64)
% Here, nf = 2*nc - 1.
% We have 64 coarse points.
% Coarse point 1 from fine points 1, 2, 3 (or 1,2)
% Coarse point i from fine points 2i-1, 2i, 2i+1.
% The range of fine points is [1, nf].

% Revised logic for general restriction operators based on common practices:

% For a "standard" multigrid coarsening (e.g., cell-centered or nodal)
% The number of coarse grid points is usually (N_fine - 1)/2 + 1 if N_fine is odd,
% or N_fine/2 if N_fine is even.

% Based on your N = [127, 64, 32], it implies:
% 127 -> 64: (127+1)/2 = 64
% 64 -> 32: 64/2 = 32

% This suggests two different restriction approaches depending on the fine grid size:
% 1. If fine grid size is odd (127): Full-weighting (1/4, 1/2, 1/4 or 1/2, 1, 1/2).
%    For '1/2, 1, 1/2' type, it usually implies that points at 2i-1 and 2i+1 are
%    boundary-like or handled specially.
%    If the stencil is applied as in your code (1/2, 1, 1/2), then for coarse point `i`:
%    `u_coarse(i) = 0.5 * u_fine(2*i-1) + 1.0 * u_fine(2*i) + 0.5 * u_fine(2*i+1)`
%    This requires `2*i-1` and `2*i+1` to be valid fine grid indices.

% Let's assume your (1/2, 1, 1/2) stencil is intended for a full-weighting-like operator
% where the boundary handling makes sense for your specific problem.
% The issue is with the `(:,2:end-1)` part.

if mod(nf, 2) ~= 0 % If fine grid has an odd number of points (e.g., 127)
    % This is likely for a "full weighting" or specific odd-point restriction.
    % The number of coarse points is (nf+1)/2. (e.g., (127+1)/2 = 64)
    % Initialize I_h_2h to be nc x nf (64 x 127)
    % For a standard 1D full weighting, the weights are often (1/4, 1/2, 1/4)
    % Your code uses (1/2, 1, 1/2) as the non-normalized weights. Let's stick to that.
    % To apply (1/2, 1, 1/2) without index issues, you need to be careful.
    % The simplest way is to ensure you don't go out of bounds.

    for i = 1:nc % Iterate through coarse grid points
        % The corresponding fine grid point is 2*i.
        % The stencil applies to 2*i-1, 2*i, 2*i+1.

        % Handle boundary conditions for the fine grid indices
        % Left boundary: only 2*i and 2*i+1 contribute if 2*i-1 is out of bounds (i=1)
        if (2*i - 1) >= 1
            I_h_2h(i, 2*i - 1) = 0.5 * c; % Apply custom 'c' factor
        end
        % Middle point (always valid if i is valid)
        if (2*i) <= nf
            I_h_2h(i, 2*i) = 1.0;
        end
        % Right boundary: only 2*i-1 and 2*i contribute if 2*i+1 is out of bounds
        if (2*i + 1) <= nf
            I_h_2h(i, 2*i + 1) = 0.5 * c;
        end
    end
    IhH = I_h_2h; % Do NOT trim columns here. The matrix should be nc x nf.

else % If fine grid has an even number of points (e.g., 64)
    % This is likely for an "injection" or "half-weighting" restriction.
    % The number of coarse points is nf/2. (e.g., 64/2 = 32)
    % Initialize I_h_2h to be nc x nf (32 x 64)

    % Your original implementation for this case looks like
    % a simple injection or a 2-point average.
    % I_h_2h(i,2*i-1)=1;
    % I_h_2h(i,2*i)=1;
    % This would mean u_coarse(i) = u_fine(2i-1) + u_fine(2i)
    % If it's half-weighting, it should be 0.5 for each.
    % Assuming your intention was this 2-point averaging with '1's.
    for i = 1:nc
        I_h_2h(i, 2*i - 1) = 1;
        I_h_2h(i, 2*i) = 1;
    end
    IhH = I_h_2h; % Already correct for this case.
end