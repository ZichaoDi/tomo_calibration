function [AFun, AtFun] = getAFunAtFun(S0_rec)
% A is diagonal matrix for tomography distortion. Each scan output does NOT depend on other scans locations
%   S:  the Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
%         AFun;         % Forward function handle:  AFun(f)  = reshape(A*f(:), gSz)
%         AtFun;        % 'Backward' function handle:  AtFun(g) = reshape(A'*g(:), fSz), where A' is conjugate transpose and A'=A.' (non-conjugate transpose) for real matrix.

[NTheta, NTau] = size(S0_rec);
S_sz = [NTheta, NTau];

SDiffCtr = zeros(size(S0_rec),'like',S0_rec); n = 2:(NTau-1); SDiffCtr(:,n) = (S0_rec(:,n+1) - S0_rec(:,n-1))/2; % SDiffCtr = [zeros(NTheta, 1), (S(:, 3:end) - S(:, 1:end-2))/2, zeros(NTheta, 1)];

AFun = @(D) reshape(SDiffCtr(:).*D(:), S_sz); % AFun  = @(D) reshape(A*D(:), S_sz); % where D is 2D drift matrix

AtFun = AFun;
end

