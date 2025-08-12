clear, clc
rng(0)

load('house.mat')

n = size(X, 1);
N = n^2;

% % first generate the blurring operator
optbl = PRblurshake('defaults');
% optbl = PRset(optbl, 'trueImage', X, 'CommitCrime', 'on');
optbl = PRset(optbl, 'trueImage', X,'CommitCrime', 'on', 'BlurLevel', 'mild');
[A1, b1, x, ProbInfo] = PRblurshake(optbl);
% A1 = speye(N);
% then generate the undersampling operator
%Undersampling ratioI
p = 0.6; % p = 0.4; % p = 0.3; % p=0.5;
%Mask
% mask = rand(n) < p; 
% ind = find(mask==1);
% M = numel(ind);
% %Mask matrix (sparse matrix in matlab)
% S = sparse(1:M, ind, ones(M, 1), M, N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we want more stuctured patches %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masktemp = rand(n);
mask = zeros(n); mask(masktemp<p) = 1;
window = zeros(21); hsw = 2; window(11-hsw:11+hsw, 11-hsw:11+hsw) = 1/(hsw+1)^2;
mask1 = conv2(mask,window); hmask = 1.7; mask1 = mask1(11:(size(mask1,1)-10),11:(size(mask1,2)-10)); mask1(mask1<=hmask) = 0; mask1(mask1>hmask) = 1;
ind = find(mask1==1);
M = numel(ind);
%Mask matrix (sparse matrix in matlab)
S = sparse(1:M, ind, ones(M, 1), M, N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finally generate the blur-inpaint operator

A = @(xx,tflag) OPblurinpaint(xx,A1,S,tflag);

% generate data
b = A_times_vec(A, x);
bdisp = zeros(N,1); bdisp(ind) = b;
undersampling = nnz(bdisp)/N;
nl = 1e-2; % nl = 0; % 
bn = PRnoise(b, nl);

reg = 0;

opth.RegParam = reg;
opth.NoiseLevel = nl;
opth.NoStop = 'on';
opth.x_true = x;
opth.DecompOut = 'on';

thr = 10^-3;

%%

% % standard LSQR
nit = 100;
[X_lsqr, info_lsqr] = IRhybrid_lsqr(A, bn, 1:nit, opth);

%%
% LR-FLSQR 
rktrunc = 20;
k = 100; eta=1.01; x0 = zeros(N,1);
parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
parameters.thr = 10^-1;
[X_flsqr, Enrm_flsqr, Rnrm_flsqr, rkX_flsqr, lambdahis_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',n,'iftrun',1,'tau',0,'reg',0,'eta',eta,'nl',nl);
% [X_flsqrvar, Enrm_flsqrvar, Rnrm_flsqrvar, rkX_flsqrvar] = LRlsqr(A, bn, x0, x, parameters);

%%
% rktrunc = 10;
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr2, Enrm_flsqr2, Rnrm_flsqr2, rkX_flsqr2] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 5;
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr3, Enrm_flsqr3, Rnrm_flsqr3, rkX_flsqr3] = LRlsqr(A, bn, x0, x, parameters);
% 
% figure, semilogy(Enrm_flsqr), hold on, semilogy(Enrm_flsqr2), semilogy(Enrm_flsqr3)
%%
% FLSQR-NNR
optnnr.p = 1;
optnnr.maxIt = 100;
optnnr.reg = reg;
optnnr.eta = eta;
optnnr.nl = nl;
optnnr.svdbasis = 'v';
optnnr.gamma = 10^-10;
optnnr.thr = 10^-3;
[X_flsqrnnr1, Enrm_flsqrnnr1, Rnrm_flsqrnnr1] = flsqr_nnr(A, bn, x0, x, optnnr);
%%
optnnr.svdbasis = 'x';
[X_flsqrnnr5, Enrm_flsqrnnr5, Rnrm_flsqrnnr5] = flsqr_nnr(A, bn, x0, x, optnnr);

%%
optnnr.maxIt = 25;
optnnr.cycles = 4;
optnnr.thr = 1e-3;
optnnr.weigthtype = 'sqrt'; % parameters.weigthtype = 'trunc'; % 
optnnr.thrstop = 1e-8;
%optnnr.reg = 0;
[X_irnlsqr_tot, X_irnlsqr_finals, X_irnlsqr_best, Enrm_irnlsqr, Rnrm_irnlsqr, ...
    cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_lsqr_nnr(A, bn, x, x0, optnnr);
%%
% optnnr.weighttype = 'trunc';
% [X_irnlsqr_finals_t, X_irnlsqr_best_t, Enrm_irnlsqr_t, Rnrm_irnlsqr_t, ...
%     cyclesIt_t, sXcycles_t, cyclesBest_t, Lambda_tot_t] = irn_lsqr_nnr(A, bn, x, x0, optnnr);


%%
%%%% displaying the results
figure
semilogy(0:length(info_lsqr.Enrm),[1;info_lsqr.Enrm],'linewidth',1.5)
hold on
semilogy(0:length(Enrm_flsqr),[1;Enrm_flsqr],'linewidth',1.5)
% semilogy(0:length(Enrm_flsqr2),[1;Enrm_flsqr2],'linewidth',1)
% semilogy(0:length(Enrm_flsqr3),[1;Enrm_flsqr3],'linewidth',1)
semilogy(0:length(Enrm_flsqrnnr1),[1;Enrm_flsqrnnr1],'linewidth',1.5)
semilogy(0:length(Enrm_flsqrnnr5),[1;Enrm_flsqrnnr5],'linewidth',1.5)
semilogy(0:length(Enrm_irnlsqr),[1;Enrm_irnlsqr],'linewidth',1.5)
legend('LSQR', 'LR-FLSQR','FLSQR-NNR(v)', 'FLSQR-NNR', 'IRN-LSQR-NNR')
legend('location','southeastoutside')
axis([0 nit 0.08 1])
xticks([0 20 40 60 80 100])
yticks([0.1 0.5 1])
set(gca,'fontsize',14)
set(gcf, 'Position',  [200, 200, 900, 300])

%%
% X_lsqr_end = reshape(info_lsqr.BestReg.X, n, n);
[~,ind] = min(info_lsqr.Enrm)
X_lsqr_end = reshape(X_lsqr(:,ind), n, n);
X_lsqr_cut = X_lsqr_end(105:105+63,124:124+63);
X_lsqr_cut = imresize(X_lsqr_cut, 4, 'nearest');
figure, imshow(X_lsqr_end, []);
figure, imshow(X_lsqr_cut, []);

[~,ind] = min(Enrm_flsqr)
X_flsqr_end = reshape(X_flsqr(:,ind), n, n);
X_flsqr_cut = X_flsqr_end(105:105+63,124:124+63);
X_flsqr_cut = imresize(X_flsqr_cut, 4, 'nearest');
figure, imshow(X_flsqr_end, []);
figure, imshow(X_flsqr_cut, []);

[~,ind] = min(Enrm_flsqrnnr1)
X_fnnr1_end = reshape(X_flsqrnnr1(:,ind), n, n);
X_fnnr1_cut = X_fnnr1_end(105:105+63,124:124+63);
X_fnnr1_cut = imresize(X_fnnr1_cut, 4, 'nearest');
figure, imshow(X_fnnr1_end, []);
figure, imshow(X_fnnr1_cut, []);

[~,ind] = min(Enrm_flsqrnnr5)
X_fnnr5_end = reshape(X_flsqrnnr5(:,ind), n, n);
X_fnnr5_cut = X_fnnr5_end(105:105+63,124:124+63);
X_fnnr5_cut = imresize(X_fnnr5_cut, 4, 'nearest');
figure, imshow(X_fnnr5_end, []);
figure, imshow(X_fnnr5_cut, []);

[~,ind] = min(Enrm_irnlsqr)
X_irnlsqr_end = reshape(X_irnlsqr_tot(:,ind), n, n);
X_irnlsqr_cut = X_irnlsqr_end(105:105+63,124:124+63);
X_irnlsqr_cut = imresize(X_irnlsqr_cut, 4, 'nearest');
figure, imshow(X_irnlsqr_end, []);
figure, imshow(X_irnlsqr_cut, []);

%%
B = reshape(bdisp, n, n);
B_cut = B(105:105+63,124:124+63);
B_cut = imresize(B_cut, 4, 'nearest');
figure, imshow(B, []);
figure, imshow(B_cut, []);

%%
X_cut = X(105:105+63,124:124+63);
X_cut = imresize(X_cut, 4, 'nearest');
figure, imshow(X, []);
figure, imshow(X_cut, []);

%%
s = svd(X); s = s/s(1);
figure, semilogy(s,'linewidth',1.5)
axis([1,256,0,1])
set(gca, 'fontsize',16)