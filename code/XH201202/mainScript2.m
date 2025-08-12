%% Chang's documentation for Wendy

% change driftType: line 55
% change maxDrift for type 1: line 59
% change noise level: line 104
% change initial lambda: line 168

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
    sampleName = 'Brain'; % choose from {'Phantom', 'Brain', ''Golosio', 'circle', 'checkboard', 'fakeRod'}; not case-sensitive
end
Nx = 256; Ny = Nx; WSz = [Ny, Nx]; %  XH: Ny = 100 -> 10 or 50 for debug.  currently assume object shape is square not general rectangle
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
NTheta = 90;  % 30/45/60/90 sample angle #. Use odd NOT even, for display purpose of sinagram of Phantom. As even NTheta will include theta of 90 degree where sinagram will be very bright as Phantom sample has verical bright line on left and right boundary.
NTau = ceil(sqrt(sum(WSz.^2))) + ceil(maxMaxDriftAll)*2; NTau = NTau + rem(NTau-Ny,2); % number of discrete beam, enough to cover object diagonal, plus tolarence with maxDrift, also use + rem(NTau-Ny,2) to make NTau the same odd/even as Ny just for perfection, so that for theta=0, we have sum(WGT, 2)' and  S(1, (1:Ny)+(NTau-Ny)/2) are the same with a scale factor
SSz = [NTheta, NTau];
driftGTAll = {0, []}; % Two Scans: 1st one no drift error. 

%% Generate drift error in pixel units for 2nd Scan
rng('default'); % same random number for repeatable test
driftType = 1;
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
%%TODO: combine types (1) + (2) + (3)
%% Generate LAll, SAll
NScan = numel(driftGTAll);
LAll = cell(2,1); 
SAll = zeros(NTheta, NTau, NScan);
if ~batchMode
    gaussianSTD = 0.01; % 0.01; % [0, 0.01, 0.05]
end
for idxScan = 1: NScan
        L = XTM_Tensor_XH(WSz, NTheta, NTau, driftGTAll{idxScan}, WGT);
        if idxScan == 1
            LNormalizer = full(max(sum(L,2))); % Compute the normalizer 0.0278 for WSz=100*100, NTheta,NTau=45*152 for 1st scan of no drift so that maximum row sum is 1 rather than a too small number
        end
        L = L/LNormalizer;
        LAll{idxScan} = L; % Matlab may only support 2D sparse matrix, must use cell not 3D array like LAll(:,:,idxScan) = L
        S = reshape(LAll{idxScan}*WGT(:), NTheta, NTau);
        % add noise to sinogram
        SMax = max(S(:));        
        rng('default'); % same random number for initial test
        SAll(:,:,idxScan) = imnoise(S/SMax, 'gaussian', 0, gaussianSTD^2)*SMax; % add noise to [0, 1 ] grayscale image, imnoise add gaussian noise but also seems to still threshold the noisy image to [0, 1]                 
%         Sdiff = SAll(:,:,idxScan) - S;
%         normeps1 = norm(Sdiff(:));
end

matfile = '20190711_phantom.mat';  % choose current name if need to save
if ~exist(matfile, 'file'); save(matfile, 'SAll', 'driftGT', 'maxDrift', 'WGT', 'LAll', 'gaussianSTD'); end

%% Show W and S
figNo = 10;
figure(figNo+1); imagesc(WGT); axis image; title('Object Ground Truth');
figure(figNo+2); multAxes(@imagesc, {SAll(:,:,2), SAll(:,:,1)}); multAxes(@title, {'Sinogram-Measured-w-Drift', sprintf('Sina-Ideal-No-Drift, psnr=%.2fdB', difference(SAll(:,:,1), SAll(:,:,2)))}); linkAxesXYZLimColorView(); 
%tilefigs;

%%
%% II: Solve linear inverse problem: L*W = S
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 500; % 100 is not enough
regWt = 1.3e-5; % 1.3e05 or 1e-5 similar or 1e-8*max(WGT(:))/LNormalizer^2;  1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxSVDofA = 1e-4; % for TV, maxSVDofA = 1e-5/1e-4 or 1e-6/LNormalizer^2 NOT maxSVDofA = 1e-6/LNormalizer,  make sure it is actually much smaller than svds(L, 1), so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};


%% Naive Reconstruction (without our joint calibration & reconstruction) using solver from XH, with L1/TV regularizer.
% Need 100/500/1000+ iteration to get good/very good/excellent result with small regularizer.
% choose small maxSVDofA to make sure initial step size is not too small. 1.8e-6 and 1e-6 could make big difference for n=2 case

solveNaive4ThreeCases = true;
if solveNaive4ThreeCases
    % solve for three cases. 
    cases = {'no drift', 'known drift', 'unknown drift, used forward model of no drift'};
    idxS = [1, 2, 2]; % corresponding S matrix for three cases
    idxL = [1, 2, 1]; % corresponding L matrix for three cases, note in the last case, the drift is unkown, we use model of zero drift

    WUncalib = zeros(WSz(1), WSz(2), numel(cases), 'like', WGT);
    for n = 1:numel(cases)      
        WUncalib(:,:,n) = solveTwist(SAll(:,:,idxS(n)), LAll{idxL(n)}, paraTwist{:});
    end
    figNo = 50;
    for n = 1:numel(cases)    
        %figure(figNo+n); W = WUncalib(:,:,n);  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWt, maxIterA));
        figure(figNo+n); W = WUncalib(:,:,n); multAxes(@imagesc, {WGT, W}); multAxes(@title, {'Ground Truth', sprintf('Rec with %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWt, maxIterA)}); linkAxesXYZLimColorView();
        PSNR_baseline = psnr(W,WGT);
        SSIM_baseline = ssim(W,WGT);
    end
    tilefigs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Our Main Algorithm to Calibrate Individual Type I/II/III drift and reconstruct object.  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sm = SAll(:,:,2);  driftGT = driftGTAll{2}; 
normeps = norm(LAll{2}*WGT(:)-Sm(:));
paraTwist{6} = 1.3e-3; % initial lambda_0
paraTwist{16} = 1e-4;
if driftType == 1
    [WRecs, SRecs, L, drift, driftAll, info] = recTompographyWithDrift(Sm, paraTwist, [], maxDrift, driftGT,'new',normeps);
    [WRecs_old, SRecs_old, L_old, drift_old, driftAll_old, info_old] = recTompographyWithDrift(Sm, paraTwist, [], maxDrift, driftGT,'old',normeps);
elseif driftType == 2    
    [WRec, S_shiftback, drift, driftAll] = recTompographyWithDrift_rotation(Sm, paraTwist, WSz, driftGT);
elseif driftType == 3
    [WRecs, SRecs, drift, driftAll,info] = recTompographyWithDrift_TypeIII(Sm, paraTwist, WSz, maxDrift, driftGT,'new',normeps);
    [WRecs_old, SRecs_old, drift_old, driftAll_old,info_old] = recTompographyWithDrift_TypeIII(Sm, paraTwist, WSz, maxDrift, driftGT,'old',normeps);
elseif driftType == 4
    error('not implemented yet!')% loop of type I and type II
end


%% create convergence plots
figure, plot(1:length(info.misfit_history),info.misfit_history,'linewidth',1.5)
hold on 
plot(1:length(info.misfit_history),info_old.misfit_history,'linewidth',1.5)
xlim([1,10])
legend('new','old')
title('misfit term history')
set(gca,'fontsize',16)   


figure, plot(1:length(info.psnr_history),info.psnr_history,'linewidth',1.5)
hold on, plot(1:length(info.psnr_history),info_old.psnr_history,'linewidth',1.5)
xlim([1,10])
legend('new', 'old')
title('PSNR history')
set(gca,'fontsize',16)   

figure, plot(1:length(info.lambda_history),info.lambda_history,'linewidth',1.5)
hold on 
plot(1:length(info.lambda_history),info_old.lambda_history,'linewidth',1.5)
xlim([1,10])
legend('new','old')
title('lambda history')
set(gca,'fontsize',16)   


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Extra Analysis for Different Interpolation Methods %%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isAnalysis = false;
% if isAnalysis
%     %% Forward Error Analysis: how good does the sinogram computed from interpolation given drift, e.g., compare P*L*w with PDriftedGT*w.   TODO: merge the test function in interpSino()
%     driftGT = driftGTAll{2};
%     SInterp1 = interpSino(SAll(:,:,1), driftGT, 'linear');   psnr1 = difference(SAll(:,:,2), SInterp1);
%     SInterp2 = interpSino(SAll(:,:,1), driftGT, 'quadratic'); psnr2 = difference(SAll(:,:,2), SInterp2); 
%     figNo = 30;
%     figure(figNo+1); multAxes(@imagesc, SAll(:,:,2), SInterp1, SInterp2); multAxes(@title, 'Drifted-GT', sprintf('Drifted-Linear-Interp,psnr=%.2fdB', psnr1), sprintf('Drifted-Quad-Interp,psnr=%.2fdB', psnr2));
% 
%     %% Backward Error Analysis: How good is calclating drift given the true sinogram pairs. TODO: merge the test function in calcDrift()
%     driftGT = driftGTAll{2};
%     interpMethods = {'linear', 'quadratic'}; 
% 
%     driftComputed = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, interpMethods{1});
%     % zero out the ground truth for the begin and end empty beams. 
%     iStart = find(cumsum(driftComputed), 1, 'first');
%     iEnd = NTau + 1 - find(cumsum(driftComputed(end:-1:1)), 1, 'first');
%     idxHasDrift = iStart:iEnd;
%     idxNoDrift = [1:iStart-1, iEnd+1:NTau];
% 
%     driftGT2 = driftGT;         driftGT2(idxNoDrift) = 0;
%     figure(1001), clf; plot([driftGT2, driftComputed], '-x'); title(sprintf('psnr=%.2fdB', difference(driftComputed, driftGT2))); legend('GT', 'rec-testInterp'); grid on; grid minor; %driftComputed-driftGT2, 'dif'
% 
% 
%     %% How good is our interpolcation model, 
%     % even given exact drift, our interpolation model is not accurate enough to produce very good reconstruction. PSNR = 24/25dB 
%     % while use the exact model, PSNR = 34dB/50dB+ use Wendy or TwIST
%     driftGT = driftGTAll{2};
%     PGT =  genDriftMatrix(driftGT, NTheta, NTau);
%     L =  PGT*LAll{1};
%     %difference(L, LAll{2})  % psnr/psnrSparse = 30.16dB/9.63dB
%     difference(L*WGT(:), LAll{2}*WGT(:)) % psnr/psnrSparse = 34.82dB/32.31dB for phantom, but 45.34dB/44.80dB for brain which is more smooth
%     figNo = 70;
%     figure(figNo+1); W2 = solveTwist(SAll(:,:,2), L, paraTwist{:}); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Interpolated Forward Model', difference(W, WGT), regType, regWt, maxIterA));
%     figure(figNo+2); W2 = solveWendy(SAll(:,:,2), L, WGT); W = W2;  imagesc(W); axis image; title(sprintf('Rec Wendy %s, PSNR=%.2fdB', 'Use Interpolated Forward Model', difference(W, WGT)));
%     figure(figNo+3); W2 = solveTwist(SAll(:,:,2), LAll{2}, paraTwist{:}); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Exact Forward Model from Drift', difference(W, WGT), regType, regWt, maxIterA));
% end
% 
% %% Our Main Algorithm to Calibrate and reconstruct object. (Type I drift)
% % Try different regWt (and use one that is in favor of other people's method) or just setup one
% isSearch4BestRegWt = false; % find best regWt for non-drift groundtruth, and use it for all
% if false
%     if isSearch4BestRegWt
%         regWtAllUncalib = 1e-5*10.^(-3:3);
%         NWt = numel(regWtAllUncalib);
%         WUncalibAllWt = zeros(WSz(1), WSz(2), numel(cases), NWt, 'like', WGT);
%         psnrAllWtUncalib = zeros(numel(cases), NWt);
%         WUncalibBestWt = zeros(WSz(1), WSz(2), numel(cases), 'like', WGT);
%         regWtBestUncalib = zeros(numel(cases), 1);
%         for n = 1:numel(cases)
%             for nWt = 1:NWt        
%                 paraTwistTry = {'xSz', WSz, 'regFun', regType, 'regWt', regWtAllUncalib(nWt), 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};
%                 WUncalibAllWt(:,:,n,nWt) = solveTwist(SAll(:,:,idxS(n)), LAll{idxL(n)}, paraTwistTry{:});
%                 psnrAllWtUncalib(n, nWt) = psnr(WUncalibAllWt(:,:,n,nWt), WGT, 1);  % WGT is normalized with maximum as 1
%             end
%             [~, idx] = max(psnrAllWtUncalib(n,:));
%             regWtBestUncalib(n) = regWtAllUncalib(idx);
%             WUncalibBestWt(:,:,n) = WUncalibAllWt(:,:,n,idx);    
%         end
%         %regWtBest = max(psnrAll
%         figNo = 100;
%         for n = 1:numel(cases)    
%             figure(figNo+n); W = WUncalibBestWt(:,:,n);  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWtBestUncalib(n), maxIterA));
%         end
%         tilefigs;
%         % setup the regWt for our calibrated algorithm
%         regWt = regWtBestUncalib(3); %phantom 100*100, maxDrift=1, regWtBestUncalib=[1e-6, 1e-6, 1e-4] % use regWtBest(2) where drift is known, or regWtBest(1) where drift is 0 
%         idx = find(strcmp(paraTwist, 'regWt'))+1;  paraTwist{idx} = regWt; % update regWt
%     end
% 
%         % Find the drift and Reconstruct by Alternative Minimization Method
%     Sm = SAll(:,:,2);  driftGT = driftGTAll{2}; 
%     [WRec, SRec, L, drift, driftAll] = recTompographyWithDrift(Sm, paraTwist, [], maxDrift, driftGT);
% end
% 
% %%
% if ~batchMode && false
% %driftIter = driftOld;
%     save('ret_20181121.mat', 'WRec', 'driftAll', 'driftGT', 'driftIter', 'maxDrift', 'SAll');
%     %%
%     paper_figures
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Recycle Code below %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Reconstruction with solver from Wendy, just least square, no regularizer, not useful
%{
WWendy = zeros(Ny, Nx, numel(cases), 'like', WGT);
for n = 1:numel(cases)
    WWendy(:,:,n) = solveWendy(SAll(:,:,idxS(n)), LAll{idxL(n)}, WGT);
end
figNo = 30;
for n = 1:numel(cases)    
    figure(figNo+n); W = WWendy(:,:,n);  imagesc(W); axis image; title(sprintf('Rec Wendy %s, PSNR=%.2fdB', cases{n}, difference(W, WGT)));
end
tilefigs;
%}

%% Wendy-L2 regularizer. This is not very useful
%{
n = 3;
WWendyL2 = solveWendyL2Reg(SAll(:,:,idxS(n)), LAll{idxL(n)}, WGT);
figure(figNo+numel(cases)+1); W = WWendyL2;  imagesc(W); axis image; title(sprintf('Rec Wendy L2-regularizer %s, PSNR=%.2fdB', cases{n}, difference(W, WGT)));
%}

%% Developping Stage: Two steps reconstruction, first TV on sinogram, then TV on object, just 1dB improvement for distorted one, not very useful
%Reconstruction with solver from XH, with L1/TV regularizer.
%{
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
regWt = 1e-8*max(WGT(:)); % 1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxIterA = 500; % 100 is not enough
maxSVDofA = 1e-6; %svds(L, 1)*1e-4; % multiply by 1e-4 to make sure it is small enough so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA};

regWtSino = regWt*10000;
paraTwistSino = {'xSz', SSz, 'regFun', regType, 'regWt', regWtSino, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', SAll(:,:,1), 'maxSVDofA', maxSVDofA};

WUncalib = zeros(WSz(1), WSz(2), numel(cases), 'like', WGT);
SAllReg = zeros(size(SAll), 'like', SAll);


IdentityMat = speye(NTheta*NTau, NTheta*NTau);
for n = 1:size(SAll,3)
    SAllReg(:,:,n) = solveTwist(SAll(:,:,n), IdentityMat, paraTwistSino{:});  % may need to change TV just along one direction
end
for n = 1:numel(cases)
%for n = 1:1    
    WUncalib(:,:,n) = solveTwist(SAllReg(:,:,idxS(n)), LAll{idxL(n)}, paraTwist{:});
end
figNo = 60;
for n = 1:numel(cases)    
    figure(figNo+n); W = WUncalib(:,:,n);  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWt, maxIterA));
end
figNo = figNo+numel(cases);
sinoCases = {'no drift', 'drift'};
for n = 1:size(SAllReg,3)
    figure(figNo+n); SReg = SAllReg(:,:,n);  imagesc(SReg); axis image; title(sprintf('Sino %s smoothed, vs Sino No Drift PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', sinoCases{n}, difference(SReg, SAll(:,:,1)), regType, regWtSino, maxIterA));
end

tilefigs([], [], [], 3); %tilefigs(handles,resize,nRows,nCols,...)
%}

