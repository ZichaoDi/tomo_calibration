clear, clc
rng(0)

load('house')

n = size(X, 1);
N = n^2;

% % first generate the blurring operator
optbl = PRblurshake('defaults');
% optbl = PRset(optbl, 'trueImage', X, 'CommitCrime', 'on');
optbl = PRset(optbl, 'trueImage', X, 'BC','reflective','CommitCrime', 'on', 'BlurLevel', 'mild');
[A1, b1, x, ProbInfo] = PRblurshake(optbl);
% A1 = speye(N);
% then generate the undersampling operator
%Undersampling ratio
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

thr = 10^-3;
k = 100; eta=1.01; x0 = zeros(N,1);

%% standard lsqr
% k = 100; eta=1.01; x0 = zeros(N,1);
% regspace1_lsqr = logspace(-10,5,16)';
% minarr1_lsqr = zeros(length(regspace1_lsqr),1);
% 
% for i = 1:length(regspace1_lsqr)
%     reg = regspace1_lsqr(i);
%     parameters = struct('maxIt',k','kappaB',1,'kappa',1,'iftrun',0,'tau',0,'reg',reg,'eta',eta,'nl',nl);
%     [X_lsqr, Enrm_lsqr, Rnrm_lsqr, rkX_lsqr, lambdahis_lsqr] = LRlsqr(A, bn, x0, x, parameters);
%     minarr1_lsqr(i) = min(Enrm_lsqr);
% end

%%
% regspace2_lsqr = linspace(0,0.002,21)';
% minarr2_lsqr = zeros(length(regspace2_lsqr),1);
% 
% for i = 1:length(regspace2_lsqr)
%     reg = regspace2_lsqr(i);
%     parameters = struct('maxIt',k','kappaB',1,'kappa',1,'iftrun',0,'tau',0,'reg',reg,'eta',eta,'nl',nl);
%     [X_lsqr, Enrm_lsqr, Rnrm_lsqr, rkX_lsqr, lambdahis_lsqr] = LRlsqr(A, bn, x0, x, parameters);
%     minarr2_lsqr(i) = min(Enrm_lsqr);
% end

%% no regularization
parameters = struct('maxIt',k','kappaB',1,'kappa',1,'iftrun',0,'tau',0,'reg',0,'eta',eta,'nl',nl);
[X_lsqr, Enrm_lsqr, Rnrm_lsqr, rkX_lsqr, lambdahis_lsqr] = LRlsqr(A, bn, x0, x, parameters);

%% 0.0011 (best) for lsqr
parameters = struct('maxIt',k','kappaB',1,'kappa',1,'iftrun',0,'tau',0,'reg',0.0011,'eta',eta,'nl',nl);
[X_lsqrreg, Enrm_lsqrreg, Rnrm_lsqrreg, rkX_lsqrreg, lambdahis_lsqrreg] = LRlsqr(A, bn, x0, x, parameters);

%% test discrep using IRTools code -- same results
% opth.RegParam = 'discrep';
% opth.NoiseLevel = nl;
% % opth.MaxIter = k;
% opth.NoStop = 'on';
% opth.x_true = x;
% opth.DecompOut = 'on';
% 
% % standard LSQR
% nit = 100;
% [X_lsqrdiscrep, info_lsqrdiscrep] = IRhybrid_lsqr(A, bn, 1:nit, opth);
%% discrep
parameters = struct('maxIt',k','kappaB',1,'kappa',1,'iftrun',0,'tau',0,'reg','discrep','eta',eta,'nl',nl);
[X_lsqrdiscrep, Enrm_lsqrdiscrep, Rnrm_lsqrdiscrep, rkX_lsqrdiscrep, lambdahis_lsqrdiscrep] = LRlsqr(A, bn, x0, x, parameters);

%% optimal
parameters = struct('maxIt',k','kappaB',1,'kappa',1,'iftrun',0,'tau',0,'reg','opt','eta',eta,'nl',nl);
[X_lsqropt, Enrm_lsqropt, Rnrm_lsqropt, rkX_lsqropt, lambdahis_lsqropt] = LRlsqr(A, bn, x0, x, parameters);

%% make plots

% figure, semilogy(0:length(Enrm_lsqr),[1;Enrm_lsqr],'linewidth',1.5), hold on, 
% semilogy(0:length(Enrm_lsqrreg),[1;Enrm_lsqrreg],'linewidth',1.5)
% semilogy(0:length(Enrm_lsqrdiscrep),[1;Enrm_lsqrdiscrep],'linewidth',1.5)
% semilogy(0:length(Enrm_lsqropt),[1;Enrm_lsqropt],'linewidth',1.5)
% legend('0','0.0011','discrep','optimal')
% axis([0 100 0.1 1])
% yticks([0.05 0.1 0.5 1])
% set(gca,'fontsize',16)


%%
% [minrelerr,ind] = min(minarr2_lsqr);
% figure,semilogy(regspace2_lsqr,minarr2_lsqr,'linewidth',1.5),hold on
% plot(regspace2_lsqr(ind),minrelerr,'p','markersize',10,'linewidth',2)
% axis([0 0.002 0.2707 0.271])
% set(gca,'fontsize',16)
% ax = gca;
% ax.XRuler.Exponent = 0;
% xticks([0,0.0005,0.001,0.0015,0.002])
% yticks([0.2707,0.2708,0.2709,0.271])
% legend('smallest relative error vs. regularization parameter')




%% flsqr-nnr

optnnr.p = 1;
optnnr.maxIt = 100;
optnnr.eta = 1.01;
optnnr.nl = nl;
optnnr.svdbasis = 'v';
optnnr.thr = 1e-3;
optnnr.gamma = 10^-10;

%%

% regspace1_nnr = logspace(-10,5,16)';
% minarr1_nnr = zeros(length(regspace1_nnr),1);
% 
% 
% for i = 1:length(regspace1_nnr)
%     optnnr.reg = regspace1_nnr(i);
%     [X_flsqrnnr1, Enrm_flsqrnnr1, Rnrm_flsqrnnr1] = flsqr_nnr(A, bn, x0, x, optnnr);
%     minarr1_nnr(i) = min(Enrm_flsqrnnr1);
% end

%%
% regspace2_nnr = logspace(-16,-12,21)';
% minarr2_nnr = zeros(length(regspace2_nnr),1);
% 
% 
% for i = 1:length(regspace2_nnr)
%     optnnr.reg = regspace2_nnr(i);
%     [X_flsqrnnr1, Enrm_flsqrnnr1, Rnrm_flsqrnnr1] = flsqr_nnr(A, bn, x0, x, optnnr);
%     minarr2_nnr(i) = min(Enrm_flsqrnnr1);
% end

%% reg parameter equal to 0 (best) for FLSQR-NNR
% optnnr.reg = 0;
% [X_flsqrnnr1reg, Enrm_flsqrnnr1reg, Rnrm_flsqrnnr1reg] = flsqr_nnr(A, bn, x0, x, optnnr);

%% discrep
optnnr.reg = 'discrep';
[X_flsqrnnr1discrep, Enrm_flsqrnnr1discrep, Rnrm_flsqrnnr1discrep] = flsqr_nnr(A, bn, x0, x, optnnr);

%% no regularization
optnnr.reg = 0;
[X_flsqrnnr1, Enrm_flsqrnnr1, Rnrm_flsqrnnr1] = flsqr_nnr(A, bn, x0, x, optnnr);

%% make plots
% figure, semilogy(0:length(Enrm_flsqrnnr1),[1;Enrm_flsqrnnr1],'linewidth',1.5), hold on, 
% % semilogy(0:length(Enrm_flsqrnnr1reg),[1;Enrm_flsqrnnr1reg],'linewidth',1.5)
% semilogy(0:length(Enrm_flsqrnnr1discrep),[1;Enrm_flsqrnnr1discrep],'linewidth',1.5)
% legend('0','discrep')
% axis([0 100 0.1 1])
% yticks([0.05 0.1 0.5 1])
% set(gca,'fontsize',16)


%%
% [minrelerr,ind] = min(minarr2_nnr);
% figure,semilogy(regspace2_nnr,minarr2_nnr,'linewidth',1.5),hold on
% plot(regspace2_nnr(ind),minrelerr,'p','markersize',10,'linewidth',2)
% axis([0 12 0.0565 0.058])
% set(gca,'fontsize',16)
% yticks([0.056,0.057,0.058])
% legend('smallest relative error vs. regularization parameter')

%% irn-lsqr-nnr
optnnr.p = 1;
optnnr.maxIt = 25;
optnnr.cycles = 4;
optnnr.thrstop = 1e-8;
optnnr.weigthtype = 'sqrt'; 
optnnr.thr = 1e-3;
optnnr.eta = 1.01;
optnnr.gamma = 10^-10;
optnnr.nl = nl;

%%
% regspace1_irn = logspace(-10,5,16)';
% minarr1_irn = zeros(length(regspace1_irn),1);
% 
% for i = 1:length(regspace1_irn)
%     optnnr.reg = regspace1_irn(i);
%     [X_irnlsqr_tot, X_irnlsqr_finals, X_irnlsqr_best, Enrm_irnlsqr, Rnrm_irnlsqr, ...
%     cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_lsqr_nnr(A, bn, x, x0, optnnr);
%     minarr1_irn(i) = min(Enrm_irnlsqr);
% end

%%
% regspace2_irn = logspace(-16,-12,21)';
% minarr2_irn = zeros(length(regspace2_irn),1);
% 
% for i = 1:length(regspace2_irn)
%     optnnr.reg = regspace2_irn(i);
%     [X_irnlsqr_tot, X_irnlsqr_finals, X_irnlsqr_best, Enrm_irnlsqr, Rnrm_irnlsqr, ...
%     cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_lsqr_nnr(A, bn, x, x0, optnnr);
%     minarr2_irn(i) = min(Enrm_irnlsqr);
% end

%% reg parameter set to 0
% optnnr.reg = 0;
% [X_irnlsqrreg_tot, X_irnlsqrreg_finals, X_irnlsqrreg_best, Enrm_irnlsqrreg, Rnrm_irnlsqrreg, ...
%     cyclesItreg, sXcyclesreg, cyclesBestreg, Lambda_totreg] = irn_lsqr_nnr(A, bn, x, x0, optnnr);

%% discrep
optnnr.reg = 'discrep';
[X_irnlsqrdiscrep_tot, X_irnlsqrdiscrep_finals, X_irnlsqrdiscrep_best, Enrm_irnlsqrdiscrep, Rnrm_irnlsqrdiscrep, ...
    cyclesItdiscrep, sXcyclesdiscrep, cyclesBestdiscrep, Lambda_totdiscrep] = irn_lsqr_nnr(A, bn, x, x0, optnnr);

%% optimal selecting stragegy
optnnr.reg = 'opt';
[X_irnlsqropt_tot, X_irnlsqropt_finals, X_irnlsqropt_best, Enrm_irnlsqropt, Rnrm_irnlsqropt, ...
    cyclesItopt, sXcyclesopt, cyclesBestopt, Lambda_totopt] = irn_lsqr_nnr(A, bn, x, x0, optnnr);

%% no regularization 
optnnr.reg = 0;
[X_irnlsqr_tot, X_irnlsqr_finals, X_irnlsqr_best, Enrm_irnlsqr, Rnrm_irnlsqr, ...
    cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_lsqr_nnr(A, bn, x, x0, optnnr);

%% make plots
% figure, semilogy(0:length(Enrm_irnlsqr),[1;Enrm_irnlsqr],'linewidth',1.5), hold on, 
% % semilogy(0:length(Enrm_irnlsqrreg),[1;Enrm_irnlsqrreg],'linewidth',1.5)
% semilogy(0:length(Enrm_irnlsqrdiscrep),[1;Enrm_irnlsqrdiscrep],'linewidth',1.5)
% semilogy(0:length(Enrm_irnlsqropt),[1;Enrm_irnlsqropt],'linewidth',1.5)
% legend('0','discrep','optimal')
% axis([0 100 0.1 1])
% yticks([0.05 0.1 0.5 1])
% set(gca,'fontsize',16)

%%
% [minrelerr,ind] = min(minarr2);
% figure,semilogy(regspace2,minarr2,'linewidth',1.5), hold on
% plot(regspace2(ind),minrelerr,'p','markersize',10,'linewidth',2)
% axis([0 200 0.065 0.07])
% set(gca,'fontsize',16)
% yticks([0.065,0.066,0.067,0.068,0.069,0.07])
% legend('smallest relative error vs. regularization parameter')
% xlabel('valuet relative error of IRN-LSQR-NNR')