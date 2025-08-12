clear, clc
rng(22)

load('peppers.mat')
load('pattern3.mat')
Xpattern = double(Xpattern);

n = size(X, 1);
N = n^2;

% % first generate the blurring operator
optbl = PRblurshake('defaults');
% optbl = PRset(optbl, 'trueImage', X, 'CommitCrime', 'on');
optbl = PRset(optbl, 'trueImage', X,'CommitCrime', 'on', 'BlurLevel', 'mild');
[A1, b1, x, ProbInfo] = PRblurshake(optbl);
% A1 = speye(N);
% finally generate the blur-inpaint operator


A = @(xx,tflag) OPblurinpaint_pattern(xx,A1,Xpattern,tflag);

% generate data
b = A_times_vec(A, x);
nl = 1e-2; % nl = 0; % 
bn = PRnoise(b, nl);

reg = 0;

opth.RegParam = reg;
opth.NoiseLevel = nl;
opth.NoStop = 'on';
opth.x_true = x;
opth.DecompOut = 'on';



%%

% % standard LSQR
nit = 100;
[X_lsqr, info_lsqr] = IRhybrid_lsqr(A, bn, 1:nit, opth);

%%
% LR-FLSQR 
rktrunc = 50;
k = 100; eta=1.01; x0 = zeros(N,1);
parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
[X_flsqr, Enrm_flsqr, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',n,'iftrun',1,'tau',0,'reg',0,'eta',eta,'nl',nl);
% [X_flsqrvar, Enrm_flsqrvar, Rnrm_flsqrvar, rkX_flsqrvar] = LRlsqr(A, bn, x0, x, parameters);

%%
% rktrunc = 10;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr1, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 30;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr3, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 40;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr4, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);

%%
% rktrunc = 20;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr5, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 60;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr6, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 70;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr7, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 80;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k,'kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr8, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
%%
% figure,semilogy(Enrm_flsqr),hold on
% semilogy(Enrm_flsqr1)
% semilogy(Enrm_flsqr5)
% semilogy(Enrm_flsqr3)
% semilogy(Enrm_flsqr4)
% legend('50','10','20','30','40')
%%
% FLSQR-NNR
optnnr.p = 1;
optnnr.maxIt = 100;
optnnr.reg = reg;
optnnr.eta = eta;
optnnr.nl = nl;
optnnr.svdbasis = 'v';
optnnr.thr = 10^-3;
optnnr.gamma = 0;
[X_flsqrnnr1, Enrm_flsqrnnr1, Rnrm_flsqrnnr1] = flsqr_nnr(A, bn, x0, x, optnnr);
%%
optnnr.svdbasis = 'x';
optnnr.gamma = 10^-10;
[X_flsqrnnr5, Enrm_flsqrnnr5, Rnrm_flsqrnnr5] = flsqr_nnr(A, bn, x0, x, optnnr);

%%
optnnr.maxIt = 25;
optnnr.cycles = 4;
optnnr.thr = 1e-3;
optnnr.weigthtype = 'sqrt'; % parameters.weigthtype = 'trunc'; % 
optnnr.thrstop = 1e-8;
optnnr.reg = 0;
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
semilogy(0:length(Enrm_flsqrnnr1),[1;Enrm_flsqrnnr1],'linewidth',1.5)
semilogy(0:length(Enrm_flsqrnnr5),[1;Enrm_flsqrnnr5],'linewidth',1.5)
semilogy(0:length(Enrm_irnlsqr),[1;Enrm_irnlsqr],'linewidth',1.5)
% semilogy(Enrm_irnlsqrreg,'linewidth',1)
% semilogy(Enrm_irnlsqr_t(1:100))
legend('LSQR', 'LR-FLSQR','FLSQR-NNR(v)', 'FLSQR-NNR', 'IRN-LSQR-NNR')
legend('location','southeastoutside')
axis([0 100 0.05 1])
yticks([0.05 0.1 0.5 1])
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
B = reshape(bn, n, n);
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