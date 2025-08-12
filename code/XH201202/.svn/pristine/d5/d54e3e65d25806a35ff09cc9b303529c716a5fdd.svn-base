function drift = calcDriftTypeIII_TVonD(Sm, S0_rec, maxDrift, regWt, driftGT, interpMethod)
%compute the drift error from sinogram matching
% Input: 
%   S0_rec:  the Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
%   Sm: the Sinogram with beamline drifted, NTheta*NTau 2D matrix, normally from measurement, we assume Sm could be interpolated from S0_rec
%   maxDrift: the maximum possible drift amount in the unit of beamlines, integer here (or pre-convert to integer: not tested)
%   driftGT:  NOT required. If provided, solver can tell how much different is our solution to groundtruth
%   interpMethod:    interplolation options
%               'linear' (same as 'central_diff'), represent sinogram s(d;w) as linear function of drift d, 2020.09.29_done
%               'piecewise-linear';  s(d; w) as pice-wise linear function of d. TODO
%    todo: add driftGT so shows optitional intermidiate improvement
% Output: 
%   drift: Type III drift. NTheta*NTau, assume small mostly [-0.5, 0.5], TODO: add check if it is too large
%          2020.07.06_still restrict drift to [-1, 1], todo change to [-0.5, 0.5]?
%          2020.08.26_working on for drift to [-n, n] where n > 1, need to
%          change solver
% See also
%   interpSinoTypeIII(); calcDrift();


if nargin < 1; testMe(); end

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

if ~exist('regWt', 'var') || isempty(regWt)
    regWt = 1e-5; % 1e-4 for cameraman?
end
if ~exist('driftGT', 'var') || isempty(driftGT)
    driftGT = [];
end

if ~exist('interpMethod', 'var') || isempty(interpMethod)
    interpMethod = 'central_diff';
end


% TODO: handle non-interger maxDrift in future
maxDriftCeil = ceil(maxDrift);
if maxDrift ~=  maxDriftCeil
    warning('non-integer maxDrift %.5f not fully tested', maxDrift)
end


%% solve via Twist with TV regularizer: 0.5*||A*D(:) - B(:)||.^2 +  regWt*TV(D) % if A is diagonal matrix, this can be solved more efficiently?
[NTheta, NTau] = size(S0_rec);
S_sz = [NTheta, NTau];

[AFun, AtFun, B] = getAFunAtFunB(Sm, S0_rec, maxDrift);
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 500; % 100 may not be enough 
paraTwist = {'xSz', S_sz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 0, 'maxIterA', maxIterA, 'xGT', driftGT, 'tolA', 1e-8, 'lb', -maxDrift, 'ub', maxDrift};

drift = solveTwist('b', B, 'AFun', AFun, 'AtFun', AtFun, paraTwist{:});

end

function testMe()
%% 2020.03_tested
% 2020.06.03_tested again for 50x50 Phantom with sinogram 45x74, drift std_=0.5
% 2020.09.29_tested for 200x200 with maxDrift = 4
% after run mainScript.m
S0_rec = SAll(:,:,1);
mask = abs(S0_rec) > eps;


%% I. Test use drifted sinogram generated from interpSinoTypeIII(), tested it is good enough for reconstruction
S2_drifted_via_interp = interpSinoTypeIII(S0_rec, driftGT, 'central_diff');  % 'linear' (default), and 'central_diff'
%S2_drifted_via_interp = interpSinoTypeIII(S0_rec, 0.5*ones(size(driftGT)));
drift = calcDriftTypeIII_TVonD(S2_drifted_via_interp, S0_rec, maxDrift, 1e-4, driftGT); % regWt=1e-4 ~ 1e-10, much larger for camera man?
figure; showDriftDifference(drift, driftGT, mask);


%% II. Test use drifted sinogram generated from measurement (with exact forward model)
Sm = SAll(:,:,2);
drift = calcDriftTypeIII_TVonD(Sm, S0_rec, maxDrift, 1e-6, driftGT); % camera man need to be much larger
figure; showDriftDifference(drift, driftGT, mask);

%% II.2. Can we use the drift from II to really reconstruct a good object,
% probably not for lena, only PSNR = 19.43dB but known gift is 65dB
% but for cameraman, can achive PSNR = 
% Under this model assumption, it is the best drift we can get as S0_rec used in II is exact SAll(:,:,1); 
L_cur = XTM_Tensor_XH(WSz, NTheta, NTau, drift, WGT);
L_cur = L_cur/LNormalizer;
W_cur = solveTwist(SAll(:,:,2), L_cur, paraTwist{:});
W = W_cur; figure; multAxes(@imagesc, {WGT, W}); multAxes(@title, {'Ground Truth', sprintf('Rec with %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'with current best drift', difference(W, WGT), regType, regWt, maxIterA)}); linkAxesXYZLimColorView();

%% IV. Test use drifted sinogram generated from exact forward model with super-resolution drift grid
N_res = 10; %2, 5, 10
%d_grid = (0:N_res-1)/N_res;
%drift_resAll = bsxfun(@plus, reshape(d_grid, N_res, 1, 1), zeros(1, NTheta, NTau));
d_grid = (-(N_res-1):N_res-1)/N_res;
N_grid = numel(d_grid);
drift_resAll = bsxfun(@plus, zeros(NTheta, NTau, 1), reshape(d_grid, 1, 1, N_grid));

L_resAll = cell(N_grid, 1);
for n = 1:N_grid
    L_resAll{n} = XTM_Tensor_XH(WSz, NTheta, NTau, drift_resAll(:, :, n));        
end
for n = 1:N_grid
    L_resAll{n} = L_resAll{n} / LNormalizer; % LNormalizer
end

S_resAll = zeros(NTheta, NTau, N_grid);
for n = 1:N_grid
    S_resAll(:, :, n) =  reshape(L_resAll{n}*WGT(:), NTheta, NTau);
end
% hack, flip to get -1, uses more memory to simplify coding for now

Sm = SAll(:,:,2);

%%

end
