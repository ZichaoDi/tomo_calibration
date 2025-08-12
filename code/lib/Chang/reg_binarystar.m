clear, clc
rng(25)
%%% best combination, so far:
% optbl = PRset(optbl, 'trueImage', 'dot2', 'BlurLevel', 'medium');
% nl = 1e-2; % 

% % first generate the blurring operator
optbl = PRblur('defaults');
optbl = PRset(optbl, 'trueImage', 'dot2', 'BlurLevel', 'medium');
% optbl = PRset(optbl, 'trueImage', 'dot2', 'BlurLevel', 'severe');
% [A, b, x, ProbInfo] = PRblurspeckle(optbl);
[A, b, x, ProbInfo] = PRblur(optbl);

N = length(x); n = sqrt(N);
nl = 1e-3; % 
% nl = 1e-2; % 
bn = PRnoise(b, nl);

reg = 0;
gamma = 10^-10;
thr = 10^-3;

k = 200; eta=1.01; x0 = zeros(N,1);


%% standard gmres
% regspace1_gmres = logspace(-10,5,16)';
% minarr1_gmres = zeros(length(regspace1_gmres),1);
% 
% for i = 1:length(regspace1_gmres)
%     reg = regspace1_gmres(i);
%     [X_gmres, Rnrm_gmres, Enrm_gmres, rkX_gmres, lambdahis_gmres] = FGMRES_LRP(A, bn, k, x0, 1, 1, x, 1, 0, reg, eta, nl);
%     minarr1_gmres(i) = min(Enrm_gmres);
% end

%%
% regspace2_gmres = logspace(-8,-6,21)';
% minarr2_gmres = zeros(length(regspace2_gmres),1);
% 
% for i = 1:length(regspace2_gmres)
%     reg = regspace2_gmres(i);
%     [X_gmres, Rnrm_gmres, Enrm_gmres, rkX_gmres, lambdahis_gmres] = FGMRES_LRP(A, bn, k, x0, 1, 1, x, 1, 0, reg, eta, nl);
%     minarr2_gmres(i) = min(Enrm_gmres);
% end

%% no regularization
[X_gmres, Rnrm_gmres, Enrm_gmres, rkX_gmres, lambdahis_gmres] = FGMRES_LRP(A, bn, k, x0, 1, 1, x, 1, 0, 0, eta, nl);

%% fixed best for gmres
reg = 1e-07;
[X_gmresreg, Rnrm_gmresreg, Enrm_gmresreg, rkX_gmresreg, lambdahis_gmresreg] = FGMRES_LRP(A, bn, k, x0, 1, 1, x, 1, 0, reg, eta, nl);
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
[X_gmresdiscrep, Rnrm_gmresdiscrep, Enrm_gmresdiscrep, rkX_gmresdiscrep, lambdahis_gmresdiscrep] = FGMRES_LRP(A, bn, k, x0, 1, 1, x, 1, 0, 'discrep', eta, nl);

%% optimal
[X_gmresopt, Rnrm_gmresopt, Enrm_gmresopt, rkX_gmresopt, lambdahis_gmresopt] = FGMRES_LRP(A, bn, k, x0, 1, 1, x, 1, 0, 'opt', eta, nl);

%% make plots

figure, semilogy(0:length(Enrm_gmres),[1;Enrm_gmres],'linewidth',1.5), hold on, 
semilogy(0:length(Enrm_gmresreg),[1;Enrm_gmresreg],'linewidth',1.5)
semilogy(0:length(Enrm_gmresdiscrep),[1;Enrm_gmresdiscrep],'linewidth',1.5)
semilogy(0:length(Enrm_gmresopt),[1;Enrm_gmresopt],'linewidth',1.5)
% semilogy(0:length(info_gmres.Enrm),[1;info_gmres.Enrm],'linewidth',1.5)
legend('0','10^{-7}','discrep','optimal')
axis([0 200 0.18 1])
yticks([0.2 0.5 1])
set(gca,'fontsize',16)


%%
[minrelerr,ind] = min(minarr2_gmres);
figure,semilogy(regspace2_gmres,minarr2_gmres,'linewidth',1.5),hold on
plot(regspace2_gmres(ind),minrelerr,'p','markersize',10,'linewidth',2)
axis([10^-8 10^-6 0.22 0.26])
set(gca,'fontsize',16)
ax = gca;
ax.XRuler.Exponent = 0;
yticks([0.22, 0.24 0.26])
xticks([10^-7 10^-6])
legend('smallest relative error vs. regularization parameter')




%% fgmres-nnr

optnnr.p = 1;
optnnr.maxIt = 200;
optnnr.regmat = 'I';
optnnr.eta = 1.01;
optnnr.nl = nl;
optnnr.gamma = 10^-10;
optnnr.thr = 10^-3;
optnnr.svdbasis = 'x';

%%

% regspace1_nnr = logspace(-10,5,16)';
% minarr1_nnr = zeros(length(regspace1_nnr),1);
% 
% 
% for i = 1:length(regspace1_nnr)
%     optnnr.reg = regspace1_nnr(i);
%     [X_fgmresnnr, Enrm_fgmresnnr, Rnrm_fgmresnnr] = fgmres_nnr(A, bn, x, x0, optnnr);
%     minarr1_nnr(i) = min(Enrm_fgmresnnr);
% end

%%
% regspace2_nnr = logspace(-10,-8,21)';
% minarr2_nnr = zeros(length(regspace2_nnr),1);
% 
% 
% for i = 1:length(regspace2_nnr)
%     optnnr.reg = regspace2_nnr(i);
%     [X_fgmresnnr, Enrm_fgmresnnr, Rnrm_fgmresnnr] = fgmres_nnr(A, bn, x, x0, optnnr);
%     minarr2_nnr(i) = min(Enrm_fgmresnnr);
% end

%% reg parameter equal to best for FGMRES-NNR
optnnr.reg = 1.584893192461111e-09;
[X_fgmresnnrreg, Enrm_fgmresnnrreg, Rnrm_fgmresnnrreg] = fgmres_nnr(A, bn, x, x0, optnnr);

%% discrep
optnnr.reg = 'discrep';
[X_fgmresnnrdiscrep, Enrm_fgmresnnrdiscrep, Rnrm_fgmresnnrdiscrep] = fgmres_nnr(A, bn, x, x0, optnnr);

%% no regularization
optnnr.reg = 0;
[X_fgmresnnr, Enrm_fgmresnnr, Rnrm_fgmresnnr] = fgmres_nnr(A, bn, x, x0, optnnr);

%% make plots
figure, semilogy(0:length(Enrm_fgmresnnr),[1;Enrm_fgmresnnr],'linewidth',1.5), hold on, 
semilogy(0:length(Enrm_fgmresnnrreg),[1;Enrm_fgmresnnrreg],'linewidth',1.5)
semilogy(0:length(Enrm_fgmresnnrdiscrep),[1;Enrm_fgmresnnrdiscrep],'linewidth',1.5)
legend('0','1.58*10^{-9}','discrep')
axis([0 200 0.2 1])
yticks([0.2 0.5 1])
set(gca,'fontsize',16)


%%
[minrelerr,ind] = min(minarr2_nnr);
figure,semilogy(regspace2_nnr,minarr2_nnr,'linewidth',1.5),hold on
plot(regspace2_nnr(ind),minrelerr,'p','markersize',10,'linewidth',2)
axis([10^-10 10^-8 0.263 0.28])
set(gca,'fontsize',16)
ax = gca;
ax.XRuler.Exponent = 0;
xticks([10^-9 10^-8])
yticks([0.265, 0.27 0.28])
legend('smallest relative error vs. regularization parameter')

%% irn-gmres-nnr
optnnr.p = 1;
optnnr.maxIt = 50;
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
%     [X_irngmres_tot, X_irngmres_finals, X_irngmres_best, Enrm_irngmres, Rnrm_irngmres, ...
%     cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_gmres_nnr(A, bn, x, x0, optnnr);
%     minarr1_irn(i) = min(Enrm_irngmres);
% end

%%
% regspace2_irn = logspace(-16,-12,21)';
% minarr2_irn = zeros(length(regspace2_irn),1);
% 
% for i = 1:length(regspace2_irn)
%     optnnr.reg = regspace2_irn(i);
%     [X_irngmres_tot, X_irngmres_finals, X_irngmres_best, Enrm_irngmres, Rnrm_irngmres, ...
%     cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_gmres_nnr(A, bn, x, x0, optnnr);
%     minarr2_irn(i) = min(Enrm_irngmres);
% end

%% reg parameter set to best which is 0
% optnnr.reg = 0;
% [X_irngmresreg_tot, X_irngmresreg_finals, X_irngmresreg_best, Enrm_irngmresreg, Rnrm_irngmresreg, ...
%     cyclesItreg, sXcyclesreg, cyclesBestreg, Lambda_totreg] = irn_gmres_nnr(A, bn, x, x0, optnnr);

%% discrep
optnnr.reg = 'discrep';
[X_irngmresdiscrep_tot, X_irngmresdiscrep_finals, X_irngmresdiscrep_best, Enrm_irngmresdiscrep, Rnrm_irngmresdiscrep, ...
    cyclesItdiscrep, sXcyclesdiscrep, cyclesBestdiscrep, Lambda_totdiscrep] = irn_gmres_nnr(A, bn, x, x0, optnnr);

%% optimal selecting stragegy
optnnr.reg = 'opt';
[X_irngmresopt_tot, X_irngmresopt_finals, X_irngmresopt_best, Enrm_irngmresopt, Rnrm_irngmresopt, ...
    cyclesItopt, sXcyclesopt, cyclesBestopt, Lambda_totopt] = irn_gmres_nnr(A, bn, x, x0, optnnr);

%% no regularization 
optnnr.reg = 0;
[X_irngmres_tot, X_irngmres_finals, X_irngmres_best, Enrm_irngmres, Rnrm_irngmres, ...
    cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_gmres_nnr(A, bn, x, x0, optnnr);

%% make plots
figure, semilogy(0:length(Enrm_irngmres),[1;Enrm_irngmres],'linewidth',1.5), hold on, 
% semilogy(0:length(Enrm_irngmresreg),[1;Enrm_irngmresreg],'linewidth',1.5)
semilogy(0:length(Enrm_irngmresdiscrep),[1;Enrm_irngmresdiscrep],'linewidth',1.5)
semilogy(0:length(Enrm_irngmresopt),[1;Enrm_irngmresopt],'linewidth',1.5)
legend('0','discrep','optimal')
axis([0 200 0.2 1])
yticks([0.2 0.5 1])
set(gca,'fontsize',16)
% set(gcf, 'Position',  [200, 200, 800, 300])
%%
% [minrelerr,ind] = min(minarr2_irn);
% figure,semilogy(regspace2_irn,minarr2_irn,'linewidth',1.5), hold on
% plot(regspace2_irn(ind),minrelerr,'p','markersize',10,'linewidth',2)
% axis([0 500 0.08655 0.0885])
% set(gca,'fontsize',16)
% yticks([0.087,0.088,0885])
% legend('smallest relative error vs. regularization parameter')
% % xlabel('value of regularization parameter')
% % ylabel('smallest relative error of IRN-LSQR-NNR')

%% compare reconstructions of hybrid GMRES and IRN-NNR (no regularization)
[~,ind] = min(Enrm_gmresopt);
X_gmres_hybrid = X_gmresopt(:,ind);
X_gmres_hybrid = reshape(X_gmres_hybrid, 256,256);

[~,ind] = min(Enrm_irngmres);
X_irngmres = X_irngmres_tot(:,ind);
X_irngmres = reshape(X_irngmres, 256,256);

band = 25; 
figure, mesh(X_gmres_hybrid(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])

figure, mesh(X_irngmres(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])

