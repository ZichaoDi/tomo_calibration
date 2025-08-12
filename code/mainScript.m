% Batch mode is run to test on different maxDrift, sampleName, and noise at 
% the same time. 
% To test on various maxDrift values, define the set "maxDriftAll" whose
% maximum is maxMaxDriftAll

%% Generate Model: 
clear

if ~exist('batchMode', 'var') || isempty(batchMode)
    batchMode = false; 
end

if ~batchMode 
    sampleName = 'Phantom';
    maxDrift = 1; % maximum drift error relative to constant scan step size
    maxMaxDriftAll = maxDrift;
    noiseLevel = 0.01;
end

% specify size of test image
Nx = 100; Ny = Nx; WSz = [Ny,Nx];

% define size of sinogram
NTheta = 45; 
NTau = ceil(sqrt(sum(WSz.^2))) + ceil(maxMaxDriftAll)*2; 
NTau = NTau + rem(NTau-Ny,2);

% generate test image, forward operator and sinogram WITHOUT drift error
driftGT0 = 0; 

% create smooth object using Gaussian filter
% WGT = createObject(sampleName,WSz);
% WGT = smooth2a(WGT,2,2);
% [~,~,x] = PRtomo(100,struct('phantomImage','smooth'));
% WGT = reshape(x,Ny,Nx);

[S0,L0,WGT,LNormalizer] = genModel(sampleName,[],WSz,NTheta,NTau,driftGT0,0);

% [L0, varargout] = build_weight_matrix_area(WGT, 0:4:176,1);
% S0 = L0*WGT(:); S0 = reshape(S0,NTheta,NTau);

% add noise
rng('default');
SMax = max(S0(:)); 
S01 = imnoise(S0/SMax, 'gaussian', 0, 0.01^2)*SMax; % noiseLevel = 0.01
S02 = imnoise(S0/SMax, 'gaussian', 0, 0.02^2)*SMax; % noiseLevel = 0.02
S03 = imnoise(S0/SMax, 'gaussian', 0, 0.03^2)*SMax; % noiseLevel = 0.03
S05 = imnoise(S0/SMax, 'gaussian', 0, 0.05^2)*SMax; % noiseLevel = 0.05

% using the same test image WGT, generate forward operator and 
% sinogram WITH drift error
rng('default');
driftGT1 = (2*rand(NTau, 1)- 1) * maxDrift;
driftGT1(1:ceil(maxMaxDriftAll)) = 0; 
driftGT1(end-ceil(maxMaxDriftAll)+1:end) = 0;
[S10,L1] = genModel('',WGT,WSz,NTheta,NTau,driftGT1,LNormalizer);

% add noise
rng('default');
SMax = max(S10(:)); 
S11 = imnoise(S10/SMax, 'gaussian', 0, 0.01^2)*SMax; % noiseLevel = 0.01
S12 = imnoise(S10/SMax, 'gaussian', 0, 0.02^2)*SMax; % noiseLevel = 0.02
S13 = imnoise(S10/SMax, 'gaussian', 0, 0.03^2)*SMax; % noiseLevel = 0.03
S15 = imnoise(S10/SMax, 'gaussian', 0, 0.05^2)*SMax; % noiseLevel = 0.05
S110 = imnoise(S10/SMax, 'gaussian', 0, 0.1^2)*SMax; % noiseLevel = 0.1


%% Reconstruction with solver from XH, with L1/TV regularizer.
% Need 100/500/1000+ iteration to get good/very good/excellent result with small regularizer.
% choose small maxSVDofA to make sure initial step size is not too small. 1.8e-6 and 1e-6 could make big difference for n=2 case
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 100; % 100 is not enough
% regWt = 1e-16;
regWtspace = logspace(-5,-1,20); % 1e-5 or 1e-8*max(WGT(:))/LNormalizer^2;  1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxSVDofA = 1e-4; % for TV, maxSVDofA = 1e-5/1e-4 or 1e-6/LNormalizer^2 NOT maxSVDofA = 1e-6/LNormalizer,  make sure it is actually much smaller than svds(L, 1), so that first step in TwIST is long enough 

res_vec = zeros(length(regWtspace),1);
wnorm_vec = zeros(length(regWtspace),1);

regWt = 0.001;
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-4};
WUncalib = solveTwist(S13, L0, paraTwist{:});

figure, imagesc(WUncalib); axis image; axis off
% title(sprintf('Baseline, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', difference(WUncalib, WGT), regType, regWt, maxIterA));


%% Reconstrucion with calibration OLD
solver = 'TwIST';
SS = S13;
noiseLevel = 0.01; optnnr.nl = norm(L1*WGT(:)-SS(:))/norm(SS(:));
regType = 'TV';
maxIter = 10;
maxIterA = 500;
maxSVDofA = 1e-4;
% regWt = noiseLevel/5;
% regWt = (0.1+maxDrift*0.01)*noiseLevel; 
% regWt = (0.001 + noiseLevel*0.02)*maxDrift;
% regWt = (0.001 + 0.03*noiseLevel)*(1+0.65*maxDrift^2)
regWt = (0.001 + 0.02*noiseLevel)*(1+0.15*maxDrift^2); 
regWt = 0.001;
% optnnr.reg = regWt; optnnr.RegParam0 = regWt;
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, ...
    'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, ...
    'maxSVDofA', maxSVDofA, 'tolA', 1e-4,'verbose',0};


[W_history_old,S_history_old,info_twist_old] = recTomoDrift_TwIST(WGT,L0,L1,SS,maxIter,solver,...
                            maxDrift,'old',regWt,1e-4,driftGT1,paraTwist,optnnr);
                      

% Reconstrucion with calibration NEW
[W_history,S_history,info_twist] = recTomoDrift_TwIST(WGT,L0,L1,SS,maxIter,solver,...
                            maxDrift,'new',regWt,0,driftGT1,paraTwist,optnnr);


% plot figures
figure, 
yyaxis left, plot(0:length(info_twist.misfit_history)-1,info_twist.misfit_history,'linewidth',1.5)
hold on 
plot(0:length(info_twist_old.misfit_history)-1,info_twist_old.misfit_history,'linewidth',1.5)
% plot(0:length(info_twist_twist.misfit_history)-1,info_twist_twist.misfit_history,'linewidth',1.5)

yyaxis right, plot(0:length(info_twist.misfit_history)-1,info_twist.psnr_history,'linewidth',1.5)

hold on, plot(0:length(info_twist_old.psnr_history)-1,info_twist_old.psnr_history,'linewidth',1.5)
% plot(0:length(info_twist_twist.psnr_history)-1,info_twist_twist.psnr_history,'linewidth',1.5)
legend('misfit new','misfit old','psnr new', 'psnr old')
% legend('misfit new','misfit dis', 'misfit twist','psnr new', 'psnr dis','psnr twist')
set(gca,'fontsize',16)      


figure,
plot(0:length(info_twist.lambda_history)-1,info_twist.lambda_history,'linewidth',1.5)
hold on, plot(0:length(info_twist_old.lambda_history)-1,info_twist_old.lambda_history,'linewidth',1.5)
legend('new lambda','old lambda')
set(gca,'fontsize',16)   


figure, imagesc(WUncalib); axis image; axis off
Wcalib = reshape(W_history_old(:,end),Ny,Nx);
figure,imagesc(Wcalib), axis image, axis off
Wcalib = reshape(W_history(:,end),Ny,Nx);
figure,imagesc(Wcalib), axis image, axis off
