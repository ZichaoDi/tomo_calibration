function drift = calcDriftTypeIII_TVonD(S2, S0_rec, driftGT, regWt, maxDrift, interpMethod)
%compute the drift error from sinogram matching
% Input: 
%   S2: the Sinogram with beamline drifted, NTheta*NTau 2D matrix, normally from measurement, we assume S2 could be interpolated from S0_rec
%   S0_rec:  the Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
%   maxDrift: the maximum possible drift amount in the unit of beamlines, integer here (or pre-convert to integer: not tested)
%   interpMethod:    interplolation options
%               'linear'
%               'quadratic', 
%    todo: add driftGT so shows optitional intermidiate improvement
% Output: 
%   drift: Type III drift. NTheta*NTau, assume small mostly [-0.5, 0.5], TODO: add check if it is too large
%          2020.07.06_still restrict drift to [-1, 1], todo change to [-0.5, 0.5]?
% See also
%   interpSinoTypeIII(); calcDrift();


if nargin < 1; testMe(); end

if ~exist('driftGT', 'var') || isempty(driftGT)
    driftGT = [];
end

if ~exist('regWt', 'var') || isempty(regWt)
    regWt = 1e-5; %
end

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

if ~exist('interpMethod', 'var') || isempty(interpMethod)
    interpMethod = 'central_diff';
end


%assert(maxDrift==floor(maxDrift), 'maxDrift need to be integer!, but it is %.5f. TODO: correct this', maxDrift)
maxDriftCeil = ceil(maxDrift);
if maxDrift ~=  maxDriftCeil
    warning('non-integer maxDrift %.5f not fully tested', maxDrift)
end
% todo: for non-interger maxDrift, add further checking, add "if alphaCur + k < maxDrift"  before "if objCur <= objOpt"

%% solve via Twist with TV regularizer: 0.5*||A*D(:) - B(:)||.^2 +  regWt*TV(D) % if A is diagonal matrix, this can be solved more efficiently?
[NTheta, NTau] = size(S0_rec);
S_sz = [NTheta, NTau];

[AFun, AtFun] = getAFunAtFun(S0_rec);
B = S2 - S0_rec;
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 500; % 100 may not be enough 
paraTwist = {'xSz', S_sz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 0, 'maxIterA', maxIterA, 'xGT', driftGT, 'tolA', 1e-4, 'lb', -maxDrift, 'ub', maxDrift};

drift = solveTwist('b', B, 'AFun', AFun, 'AtFun', AtFun, paraTwist{:});

end

function testMe()
%% 2020.03_tested
% 2020.06.03_tested again for 50x50 Phantom with sinogram 45x74, drift std_=0.5
% after run mainScript.m
S0_rec = SAll(:,:,1);
mask = (S0_rec>eps);


%% I. Test use drifted sinogram generated from interpSinoTypeIII()
S2_drifted_via_interp = interpSinoTypeIII(S0_rec, driftGT, 'central_diff');  % 'linear' (default), and 'central_diff'
%S2_drifted_via_interp = interpSinoTypeIII(S0_rec, 0.5*ones(size(driftGT)));
drift = calcDriftTypeIII_TVonD(S2_drifted_via_interp, S0_rec, driftGT, 1e-6); % 1e-5 ~ 1e-10
%drift = -calcDriftTypeIII_TVonD(S0_rec, S2_drifted_via_interp, -driftGT); % this seems to be slightly better, why?

driftGT_masked = driftGT; driftGT_masked(~mask) = 0; drift_masked = drift; drift_masked(~mask) = 0;
figure; multAxes(@imagesc, {driftGT_masked, drift_masked, drift_masked-driftGT_masked}); multAxes(@title, {'drift-GT', 'drift-Rec', sprintf('psnr=%.2fdB', difference(drift_masked, driftGT_masked))}); linkAxesXYZLimColorView(); multAxes(@colorbar);

% old 1D comparison
dGT=driftGT(mask); d=drift(mask); err_mean =  mean(abs(d(:)-dGT(:)));
figure, plot(dGT, '-rx'); hold on; plot(d, '--bo'); title(sprintf('mean error =%.2e',err_mean)); legend('GT', 'calc'); grid on; ylim([-1.1, 1.1]);
axes('Position',[.7 .7 .2 .2]); box on; histogram(d-dGT, 50); title('hist of error')


%% II. Test use drifted sinogram generated from exact forward model
S2_drfited_via_exactforwardmodel = SAll(:,:,2);
drift = calcDriftTypeIII_TVonD(S2_drfited_via_exactforwardmodel, S0_rec, driftGT, 2*1e-5);
%drift = calcDriftTypeIII_TVonD(S2_drfited_via_exactforwardmodel, S0_rec, 1,  'central_diff');  % 'linear' (default), and 'central_diff'

driftGT_masked = driftGT; driftGT_masked(~mask) = 0; drift_masked = drift; drift_masked(~mask) = 0;
figure; multAxes(@imagesc, {driftGT_masked, drift_masked, drift_masked-driftGT_masked}); multAxes(@title, {'drift-GT', 'drift-Rec', sprintf('psnr=%.2fdB', difference(drift_masked, driftGT_masked))}); linkAxesXYZLimColorView(); multAxes(@colorbar);

% old 1D comparison
dGT=driftGT; dGT(~mask)=0; d=drift; d(~mask)=0; err_mean =  mean(abs(d(:)-dGT(:)));
figure, plot(dGT, '-rx'); hold on; plot(d, '--bo'); title(sprintf('mean error =%.2e',err_mean)); legend('GT', 'calc'); grid on; ylim([-1.1, 1.1]);
axes('Position',[.7 .7 .2 .2]); box on; histogram(d-dGT, 50); title('hist of error')

%% II.2. Can we use the drift from II to really reconstruct a good object, probably not, only PSNR = 19.43dB but known gift is 65dB
% Under this model assumption, it is the best drift we can get as S0_rec used in II is exact SAll(:,:,1); 
L_cur = XTM_Tensor_XH(WSz, NTheta, NTau, drift, WGT);
L_cur = L_cur/LNormalizer;
W_cur = solveTwist(SAll(:,:,2), L_cur, paraTwist{:});
W = W_cur; figure; multAxes(@imagesc, {WGT, W}); multAxes(@title, {'Ground Truth', sprintf('Rec with %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'with current best drift', difference(W, WGT), regType, regWt, maxIterA)}); linkAxesXYZLimColorView();
end
