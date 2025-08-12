[S, S0_GT, driftGT] = readData();
[NTheta, NTau] = size(S); SSz = [NTheta, NTau];
%% Object domain, Nx = NTau/sqrt(2) - maxDrift*2
Nx = 50; Ny = Nx; WSz = [Ny, Nx]; % current assume object shape is square not general rectangle

L0 = XTM_Tensor_XH(WSz, NTheta, NTau);  % no drift
LNormalizer = full(max(sum(L0,2)));
L0 = L0/LNormalizer; % no drift
WGT = [];
%%
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 500; % 100 is not enough
regWt = 1e-5; % 1e-5 or 1e-8*max(WGT(:))/LNormalizer^2;  1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxSVDofA = 1e-4; % for TV, maxSVDofA = 1e-5/1e-4 or 1e-6/LNormalizer^2 NOT maxSVDofA = 1e-6/LNormalizer,  make sure it is actually much smaller than svds(L, 1), so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};


%%
isAnalysis = false;
if isAnalysis
    %
    W = solveTwist(S0_GT, L0, paraTwist{:});
    figure(31); clf; imagesc(W); colorbar; title('Reconstructed Object from none-drifted sino'); axis image;

    WUncalib = solveTwist(S, L0, paraTwist{:});
    figure(41); clf; imagesc(WUncalib); colorbar; title('Reconstructed Object from drifted sino'); axis image;
    tilefigs;

    %compare with S0_estimated and S to find drift?
    S0_estimated = reshape(L0*WUncalib(:), SSz);
    figure(51); clf; multAxes(@imshow, {S, S0_estimated, S0_GT}); tilefigs; linkAxesXYZLimColorView; multAxes(@title, {'S', 'S0\_estimated', 'S0\_GT'})
end


%%
%[WRec, S_shiftback, drift, driftAll] = recTompographyWithDrift_rotation(S, paraTwist, WSz, driftGT);
[WRec, S_shiftback, drift, driftAll] = recTompographyWithDrift_rotation(S, paraTwist, WSz);


%%
function [S, S0_GT, driftGT] = readData()
%% read measured sino: 
matFile = '/sandbox/Dropbox/Tomography/data/XRF_raw.mat';
XRF_raw = load(matFile, 'XRF_raw'); XRF_raw = XRF_raw.XRF_raw;
S = double(XRF_raw(:,:,4,86)); % distorted sinogram 48*81
S = flipud(S/max(S(:)));
figure(11);clf; imshow(S)

%% read none-shifted groundtruth sino S0
S0_raw = imread('/sandbox/Dropbox/Tomography/data/2019.10.22_ret.png'); 
S0_GT = sum(double(S0_raw), 3); S0_GT = S0_GT/max(S0_GT(:)); 
N = size(S0_GT,1)/2; 
S0_GT = imresize(S0_GT(N+1:end, :), size(S)); S0_GT = S0_GT/max(S0_GT(:));

%figure(21); clf; multAxes(@imshow, {S0_raw, S0_GT})

driftGT = shiftBackward(S0_GT, S, 10);
end