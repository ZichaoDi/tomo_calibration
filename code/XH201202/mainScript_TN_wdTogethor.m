%% 2018.08.17_change script to function and remove global variables. Due to a few reasons
% need to take care of the forward model matrix L with function of driftAll in scan anyway.
% not easy to debug or change parameters: e.g. change N in do_setup_XH script will raise error, may need to change multiple files
% 2nd run of programs after change parameters may have error (due to global variable).  clear all won't help and need to restart matlab.
% 2018.11.30 added new test examples
% 2018.12.03 correct a very annoying bug effect about sparse matrix.
%     Matlab doesn't support 3D sparse matrix, use LAll(:,:,1) = L will make LAll a full matrix (without any message), 
%     even though L is a sparse matrix! Use LAll
% Change to Separable: L*A - S ---> log(L*A./S), may need to add small positive number to avoid divide by 0
% Test big noise, large drift. whenver drift is visible. 
% Ratio: drift vs total size. 

%% I: Generate Model: 
% (0) in batchmode, such as generate multiple running for different maxDrift, sampleName, and noise for a paper  
if ~exist('batchMode', 'var') || isempty(batchMode); batchMode = false; end

% various pixel drift, various sample noise (maybe varies angle number)
% Angle: 20, 40, 90, 180.  
% (1) Generate Object
% sampleName: the type and name to specify a sample object, not case-sensitive
% Ny*Nx:    Sample object size
% WGT:      Sample object ground truth

if ~batchMode
    sampleName = 'Phantom'; % choose from {'Phantom', 'Brain', ''Golosio', 'circle', 'checkboard', 'fakeRod'}; not case-sensitive
end
Nx = 50; Ny = Nx; WSz = [Ny, Nx]; %  XH: Ny = 100 -> 10 or 50 for debug.  currently assume object shape is square not general rectangle
WGT = createObject(sampleName, WSz); WGT = WGT/max(WGT(:)); assert(all(WGT(:)>=0), 'Groundtruth object should be nonnegative');% Object without noise. always assume WGT0 is normailzed with maximum 1.

% (2) Scan Object: Create Forward Model and Sinograms
% NTheta*NTau:  Sinogram Size
% NScan:        Number of Scans, for example, NScan = 2, the first scan no drift, the 2nd scan has drift
% driftGTAll:  Drift amount for all scans of NScan*1 cell, driftGTAll{n} for nth scan
% LAll:         Forward Model for all scans of NScan*1 cell, LAll{n} for the nth scan
% SAll:         Sinogram for all scans of NTheta*NTau*NScan, SAll(:,:,n) for the nth scan

if ~batchMode
    maxDrift = 1; % maximum drift error; choose from [0.5, 1, 2, 4, 8] or [0.5, 1, 2, 3, 5] etc.
    maxMaxDriftAll = maxDrift; % in batch mode we consider different maxDrift from maxDriftAll set, where the maximum is maxMaxDriftAll
end
NTheta = 31;  % 30/60/45/90 sample angle #. Use odd NOT even, for display purpose of sinagram of Phantom. As even NTheta will include theta of 90 degree where sinagram will be very bright as Phantom sample has verical bright line on left and right boundary.
NTau = ceil(sqrt(sum(WSz.^2))) + ceil(maxMaxDriftAll)*2; NTau = NTau + rem(NTau-Ny,2); % number of discrete beam, enough to cover object diagonal, plus tolarence with maxDrift, also use + rem(NTau-Ny,2) to make NTau the same odd/even as Ny just for perfection, so that for theta=0, we have sum(WGT, 2)' and  S(1, (1:Ny)+(NTau-Ny)/2) are the same with a scale factor
SSz = [NTheta, NTau];
driftGTAll = {0, []}; % Two Scans: 1st one no drift error. 

%% Generate drift error in pixel units for 2nd Scan
rng('default'); % same random number for repeatable test
driftType = 3;
if driftType == 1
    % (1) Type I Drift (systematic beam drift). Assume same drift for all the angles for the same beam. NTau*1 or 1*NTau unknowns
    %TODO: change from NTau*1 to 1*NTau, so that could match Ntheta*NTau dimension in correct order
    maxDrift = 1; 
    driftGT = (2*rand(NTau, 1)- 1) * maxDrift; 
    driftGT(1:ceil(maxMaxDriftAll)) = 0;  driftGT(end-ceil(maxMaxDriftAll)+1:end) = 0; % remember NTau = diagonal + ceil(maxMaxDrift)*2 with zero padding of ceil(maxDrift)
elseif driftType == 2
    % (2) Type II Drift (systematic angle drift). Assume same drift for all the beams for the same angle. NTheta*1 unkowns 
    maxDrift = 1; 
    driftGT = (2*rand(NTheta, 1)- 1) * maxDrift; 
    driftGT(1:ceil(maxMaxDriftAll)) = 0;  driftGT(end-ceil(maxMaxDriftAll)+1:end) = 0; % remember NTau = diagonal + ceil(maxMaxDrift)*2 with zero padding of ceil(maxDrift)
elseif driftType == 3
    % (3) Type III Drift (random rotation error + beam drift) Assume drift for each angle and each scan beam position. NTheta*NTau unkowns
    maxDrift = 1;

    % (3.1) random Gaussian, truncated to [-maxDrift, maxDrift]
    std_ = maxDrift/2; % maxDrift/2 or maxDrift/3 . 95% within 2*std_, so 95% are within [-0.5, 0.5] or [-1, 1] for std_ = 0.25 or 0.5 
    driftGT = std_ *randn(NTheta, NTau); 
    driftGT(driftGT >= maxDrift) = maxDrift - eps; driftGT(driftGT <= -maxDrift ) = -maxDrift + eps; % threshold large    
    
    % (3.2) random Uniform     
    driftGT = (2*rand(NTheta, NTau)- 1) * maxDrift; 
    
    % (3.3) structured, also has L1-sparsity prior, we first try use a real image as drift        
    driftGT = imresize(imread('cameraman512.bmp'), [NTheta, NTau]);  % driftGT = imresize(WGT', [NTheta, NTau]);
    driftGT = (normalize01(driftGT)*2 - 1)*maxDrift; 
    %driftGT(driftGT >= 1) = 1 - eps; driftGT(driftGT <= -1 ) = -1 + eps; % threshold large, old code which is wrong when maxDrift not 1
    figure, imagesc(driftGT); axis image; colorbar; title('2D drift')
elseif driftType == 4
    % (4) Type I + Type II. NTheta + NTau
    maxDrift = 1; 
    driftGT1 = (2*rand(1, NTau)- 1) * maxDrift; 
    driftGT1(1:ceil(maxMaxDriftAll)) = 0;  driftGT1(end-ceil(maxMaxDriftAll)+1:end) = 0; % remember NTau = diagonal + ceil(maxMaxDrift)*2 with zero padding of ceil(maxDrift)

    maxDrift = 0; 
    driftGT2 = (2*rand(NTheta, 1)- 1) * maxDrift; 
    driftGT2(1:ceil(maxMaxDriftAll)) = 0;  driftGT2(end-ceil(maxMaxDriftAll)+1:end) = 0; % remember NTau = diagonal + ceil(maxMaxDrift)*2 with zero padding of ceil(maxDrift)
    
    driftGT = bsxfun(@plus, driftGT1, driftGT2);
end

driftGTAll{2} = driftGT;
%%TODO: combine three types
% (1) + (2) + (3)
%driftGTAll{2} = repmat(driftGTAll{2}(:)', NTheta, 1) + driftTypeII;

%% Generate LAll, SAll
NScan = numel(driftGTAll);
LAll = cell(2,1); 
SAll = zeros(NTheta, NTau, NScan);
if ~batchMode
    gaussianSTD = 0.0; % [0, 0.01, 0.05]
end
for idxScan = 1: NScan
        L = XTM_Tensor_XH(WSz, NTheta, NTau, driftGTAll{idxScan}, WGT);
        if idxScan == 1
            global LNormalizer
            LNormalizer = full(max(sum(L,2))); % Compute the normalizer 0.0278 for WSz=100*100, NTheta,NTau=45*152 for 1st scan of no drift so that maximum row sum is 1 rather than a too small number
        end
        L = L/LNormalizer;
        LAll{idxScan} = L; % Matlab may only support 2D sparse matrix, must use cell not 3D array like LAll(:,:,idxScan) = L
        S = reshape(LAll{idxScan}*WGT(:), NTheta, NTau);
        % add noise to sinogram
        SMax = max(S(:));        
        rng('default'); % same random number for initial test
        SAll(:,:,idxScan) = imnoise(S/SMax, 'gaussian', 0, gaussianSTD^2)*SMax; % add noise to [0, 1 ] grayscale image, imnoise add gaussian noise but also seems to still threshold the noisy image to [0, 1]                 
end


%% Show W and S
figNo = 10;
figure(figNo+1); imagesc(WGT); axis image; title('Object Ground Truth');
figure(figNo+2); multAxes(@imagesc, {SAll(:,:,2), SAll(:,:,1)}); multAxes(@title, {'Sinogram-Measured-w-Drift', sprintf('Sina-Ideal-No-Drift, psnr=%.2fdB', difference(SAll(:,:,1), SAll(:,:,2)))}); linkAxesXYZLimColorView(); 
%tilefigs;


%%  Wendy TN
S_measured = SAll(:,:,end);
%% (1) Reconstruction with solver from Wendy, just least square, no regularizer, with w, d tegother
wd = zeros(prod(WSz)+prod(SSz), 1);
global N NF maxiter bounds  W0 err0
maxiter=10;
bounds=1;
%%=====================================
low=[zeros(prod(WSz),1);-maxDrift.*ones(prod(SSz),1)];
up=[inf*ones(prod(WSz),1);maxDrift.*ones(prod(SSz),1)];
fctn = @(x)sfun_grad_wd_driftTypeIII(x, S_measured, WSz);
% %% ============================================
% W0=W0(1:prod(WSz));
% low=[zeros(prod(WSz),1)];
% up=[inf*ones(prod(WSz),1)];
% %% ============================================
NF = [0*N; 0*N; 0*N];
maxOuterIter=50;
for iter=1:maxOuterIter
    N=prod(WSz);
    if(iter==1)
        x0_w=wd(1:prod(WSz));
        L_d=LAll{1};
    else
        x0_w=x_w;
    end
    W0=x0_w;
    err0=norm(x0_w);
    global L_d
    fctn_w = @(x)sfun_radon_XH(x,S_measured,L_d,'LS',WSz);
    [x_w,f,g,ierror] = tnbc(x0_w,fctn_w,low(1:N),up(1:N)); % algo='TNbc';

    N=prod(SSz);
    if(iter==1)
        x0_d=wd(prod(WSz)+1:end);
    else
        x0_d=x_d;
    end
    W0=x0_d;
    err0=norm(x0_d);
    fctn_d = @(x)sfun_drift(x,x_w,S_measured,WSz);
    [x_d,f,g,ierror] = tnbc(x0_d,fctn_d,low(prod(WSz)+1:end),up(prod(WSz)+1:end)); % algo='TNbc';
end

