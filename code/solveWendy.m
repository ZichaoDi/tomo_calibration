function W = solveWendy(S, L, WGT, WInitialGuess)
% S:    2D sinogram
% L:    sparse forward model matrix
% WGT:  the 2D ground truth object which is only used for output the reconstruction error. It is vectorized and assigned to
%       globel variable W0 = WGT(:)
% WInitialGuess: the 2D intial guess of W, that is used as starting point for optimization. The default are all zeros.

% Temp fix.The real fix wouls be set the convergence critira |g| relative to data
% Bug of code, maximum S must be around 10,  hack to fix. 
% intensity
%scale = 10/max(S(:));
scale = 40; 
S = S*scale;
WGT = WGT*scale;



if ~exist('WInitialGuess', 'var') || isempty(WInitialGuess)
    WInitialGuess = ones(size(WGT), 'like', WGT);
end

global maxiter W0 err0 bounds NumElement N NF Joint
maxiter = 50; W0 = WGT(:); err0 = norm(W0(:)); bounds = 1; NumElement = 1; N = size(WGT, 1);  NF = [0*N; 0*N; 0*N]; Joint = -1; % 0: XRF; -1: XTM; 1: Joint inversion


fctn = @(x)sfun_radon_XH(x, S, L, 'LS');% on attenuation coefficients miu; % LS or EM
[x, f, g, ierror] = tnbc (WInitialGuess(:), fctn, zeros(N^2,1), inf*ones(N^2,1)); % algo='TNbc';
%[x, f, g, ierror] = tnbc (W0, fctn, zeros(N^2,1), inf*ones(N^2,1)); % algo='TNbc';
W = reshape(x, N, N);

clear maxiter W0 err0 bounds NumElement N NF Joint

% Temp bug fix: scale back
W = W/scale;
WGT = WGT/scale;

end