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

%%
X = reshape(x,256,256);
B = reshape(bn,256,256);
band = 25;
figure, imagesc(X), axis image, axis off
figure
mesh(X(128-band:128+band,128-band:128+band))
axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])

figure, imagesc(B), axis image, axis off
figure
mesh(B(128-band:128+band,128-band:128+band))
axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])

%%

%%%% SOLVERS BASED ON THE ARNOLDI ALGORITHM %%%%

% %%%% standard GMRES
maxIt = 200;
opth.RegParam = reg;
opth.NoiseLevel = nl;
opth.NoStop = 'on';
opth.x_true = x;
opth.DecompOut = 'on';

[X_gmres, info_gmres] = IRhybrid_gmres(A, bn, 1:maxIt, opth);

%%

%%%% irn GMRES
eta = 1.01;
x0 = zeros(N, 1);
p = 1;
maxIt = 50;

parameters.cycles = 4;
parameters.thr = thr;
parameters.reg = reg; % parameters.reg = 'discrep'; % 
parameters.p = p;
parameters.maxIt = maxIt;
parameters.eta = eta;
parameters.nl = nl;
parameters.weigthtype = 'sqrt'; % parameters.weigthtype = 'trunc'; % 
parameters.thrstop = 1e-8;
parameters.gamma = gamma;

[Xtot_irn, Xfinals_irn, Xbest_irn, Enrm_tot_irn, Rnrm_tot_irn, cyclesIt_irn, sXcycles_irn, cyclesBest_irn, Lambda_tot] = irn_gmres_nnr(A, bn, x, x0, parameters);


%%
parameters.reg = 'opt';
[Xtot_irnopt, Xfinals_irnopt, Xbest_irnopt, Enrm_tot_irnopt, Rnrm_tot_irnopt, cyclesIt_irnopt, sXcycles_irnopt, cyclesBest_irnopt, Lambda_totopt] = irn_gmres_nnr(A, bn, x, x0, parameters);

%%
parameters.reg = 'discrep';
[Xtot_irndiscrep, Xfinals_irndiscrep, Xbest_irndiscrep, Enrm_tot_irndiscrep, Rnrm_tot_irndiscrep, cyclesIt_irndiscrep, sXcycles_irndiscrep, cyclesBest_irndiscrep, Lambda_totdiscrep] = irn_gmres_nnr(A, bn, x, x0, parameters);
% reg.value = 1;
% parameters.reg = reg; % 
% parameters.reg = 'discrep'; % 
% [Xtot_irnp, Xfinals_irnp, Xbest_irnp, Enrm_tot_irnp, Rnrm_tot_irnp, cyclesIt_irnp, sXcycles_irnp, cyclesBest_irnp, Lambda_irnp] = irn_gmres_nnr(A, bn, x, x0, parameters);

% reg.value = 1;
% reg.rule = 'discrep';
% parameters.reg = reg; % 
% [Xfinals_irnp1, Xbest_irnp1, Enrm_tot_irnp1, Rnrm_tot_irnp1, cyclesIt_irnp1, sXcycles_irnp1, cyclesBest_irnp1, Lambda_irnp1] = irn_gmres_nnr(A, bn, x, x0, parameters);
% 
% reg.value = 1e8;
% reg.rule = 'discrepvar';
% parameters.reg = reg; %
% [Xfinals_irnp2, Xbest_irnp2, Enrm_tot_irnp2, Rnrm_tot_irnp2, cyclesIt_irnp2, sXcycles_irnp2, cyclesBest_irnp2, Lambda_irnp2] = irn_gmres_nnr(A, bn, x, x0, parameters);

%%
parameters.p = p;
parameters.reg = reg; % 
parameters.weigthtype = 'trunc'; % parameters.weigthtype = 'sqrt'; % 

[Xtot_irnt, Xfinals_irnt, Xbest_irnt, Enrm_tot_irnt, Rnrm_tot_irnt, cyclesIt_irnt, sXcycles_irnt, cyclesBest_irnt] = irn_gmres_nnr(A, bn, x, x0, parameters);

% parameters.reg = 'discrep'; % 
% [Xtot_irntp, Xfinals_irntp, Xbest_irntp, Enrm_tot_irntp, Rnrm_tot_irntp, cyclesIt_irntp, sXcycles_irntp, cyclesBest_irntp] = irn_gmres_nnr(A, bn, x, x0, parameters);
% 

%%
%%%% FGMRES, nnr 

optnnr.p = p;
optnnr.maxIt = 200;
optnnr.regmat = 'I';
optnnr.reg = 0;
optnnr.eta = eta;
optnnr.nl = nl;
optnnr.gamma = gamma;
optnnr.thr = thr;
% optnnr.svdbasis = 't';
% [X_fnnrt,Enrm_fnnrt,Rnrm_fnnrt] = fgmres_nnr(A, bn, x, x0, optnnr);
optnnr.svdbasis = 'x';
[X_fnnr,Enrm_fnnr,Rnrm_fnnr] = fgmres_nnr(A, bn, x, x0, optnnr);
optnnr.svdbasis = 'v';
[X_fnnrt,Enrm_fnnrt,Rnrm_fnnrt] = fgmres_nnr(A, bn, x, x0, optnnr);

%%
% optnnr.reg = 'discrep'; % 
% [X_fnnrp,Enrm_fnnrp,Rnrm_fnnrp,Lambda_fnnrp] = fgmres_nnr(A, bn, x, x0, optnnr);

% optnnr.p = 0;
% [X_fnnr0,Enrm_fnnr0,Rnrm_fnnr0] = fgmres_nnr(A, bn, x, x0, optnnr);
% 

%%

%%%% OTHER METHODS %%%%

maxIt = 200;

%%% LR-FGMRES
tau = 2;
[X_LRgm_2,Rnrm_LRgm_2,Enrm_LRgm_2] = FGMRES_LRP(A,bn,maxIt,x0,tau,tau,x,1e-4,1,reg,1.01,nl);

%%% RS-LR_GMRES
[X_RS_2, RelRes_RS_2, RelErr_RS_2] = RS_GMRES_LRP(A,bn,10,20,x0,tau,tau,x);

%%
%%% LR-FGMRES
tau = 30;
[X_LRgm,Rnrm_LRgm,Enrm_LRgm] = FGMRES_LRP(A,bn,maxIt,x0,tau,tau,x,1e-4,1,reg,1.01,nl);

%%
%%% RS-LR_GMRES
tau = 30;
[X_RS, RelRes_RS, RelErr_RS] = RS_GMRES_LRP(A,bn,5,40,x0,tau,tau,x);

%%
%%% SVT
[X_SVT,RelErr_SVT,RelRes_SVT] = SVT(A,bn,x,1,maxIt,1,2);

%%
% %%%% suitable plots...
% figure, semilogy(Enrm_tot_irn)
% hold on
% semilogy(Enrm_tot_irnt)
% semilogy(Enrm_fnnr)
% 
% figure, semilogy(Enrm_tot_irnp1)
% hold on
% semilogy(Enrm_tot_irnp2)
% semilogy(Enrm_fnnrp)
% 
[~,ind] = min(info_gmres.Enrm);
Xb_gm = reshape(X_gmres(:,ind), n, n); % Xb_gm = reshape(Xbest_irn(:,1), n, n); % 
s_gm = svd(Xb_gm);
Xb_irn = reshape(Xfinals_irn(:,4), n, n); % Xb_irn = reshape(Xbest_irn(:,2), n, n); % 
s_irn = svd(Xb_irn);
Xb_irnt = reshape(Xfinals_irnt(:,4), n, n); % Xb_irnt = reshape(Xbest_irnt(:,2), n, n); % 
s_irnt = svd(Xb_irnt);
[~, ind] = min(Enrm_LRgm);
Xb_lrfgm = reshape(X_LRgm(:,ind), n, n);
s_lrfgm = svd(Xb_lrfgm);
[~, ind] = min(RelErr_SVT);
Xb_svt = reshape(X_SVT(:,ind), n, n);
s_svt = svd(Xb_svt);
[~, ind] = min(Enrm_fnnr);
Xb_fnnr = reshape(X_fnnr(:,ind), n, n);
s_fnnr = svd(Xb_fnnr);
[~, ind] = min(RelErr_RS);
Xb_rs = reshape(X_RS(:,ind), n, n);
s_rs = svd(Xb_rs);

figure, semilogy(s_gm/s_gm(1),'o', 'linewidth', 1.5, 'markersize',10)
hold on
semilogy(s_irn/s_irn(1), 's', 'linewidth', 1.5, 'markersize',10)
% semilogy(s_irnt/s_irnt(1),'s', 'linewidth', 1.5, 'markersize',10)
semilogy(s_fnnr/s_fnnr(1),'d', 'linewidth', 1.5, 'markersize',10)
semilogy(s_lrfgm/s_lrfgm(1),'h', 'linewidth', 1.5, 'markersize',10)
semilogy(s_svt/s_svt(1),'*', 'linewidth', 1.5, 'markersize',10)
% semilogy(s_rs/s_rs(1),'*', 'linewidth', 1.5, 'markersize',10)
legend('GMRES','IRN-GMRES-NNR','FGMRES-NNR','LR-FGMRES','SVT')
legend('location','northeast')
axis([1 80 10^-3 1])
set(gca, 'fontsize',16)
xticks([1 20 40 60 80]), yticks([10^-3,10^-2,10^-1,1])




%% plot s at each out it of irn
X_irn1 = reshape(Xfinals_irn(:,1), n, n); s_irn1 = svd(X_irn1);
X_irn2 = reshape(Xfinals_irn(:,2), n, n); s_irn2 = svd(X_irn2);
X_irn3 = reshape(Xfinals_irn(:,3), n, n); s_irn3 = svd(X_irn3);
X_irn4 = reshape(Xfinals_irn(:,4), n, n); s_irn4 = svd(X_irn4);
X_irnt1 = reshape(Xfinals_irnt(:,1), n, n); s_irnt1 = svd(X_irnt1);
X_irnt2 = reshape(Xfinals_irnt(:,2), n, n); s_irnt2 = svd(X_irnt2);
X_irnt3 = reshape(Xfinals_irnt(:,3), n, n); s_irnt3 = svd(X_irnt3);
X_irnt4 = reshape(Xfinals_irnt(:,4), n, n); s_irnt4 = svd(X_irnt4);

figure, semilogy(s_irn1(1:20)/s_irn1(1), 'o', 'linewidth', 1.5, 'markersize',10)
hold on
semilogy(s_irn2(1:20)/s_irn2(1), 's', 'linewidth', 1.5, 'markersize',10)
semilogy(s_irn3(1:20)/s_irn3(1), '+', 'linewidth', 1.5, 'markersize',10)
semilogy(s_irn4(1:20)/s_irn4(1), 'd', 'linewidth', 1.5, 'markersize',10)
legend('k = 1','k = 2','k = 3','k = 4')
axis([1 20 10^-3 1])
set(gca, 'fontsize',16), xticks([5 10 15 20]), yticks([10^-3,10^-2,1e-1,1])


% figure, semilogy(s_irnt1(1:20)/s_irnt1(1), 'o', 'linewidth', 1.5, 'markersize',10)
% hold on
% semilogy(s_irnt2(1:20)/s_irnt2(1), 's', 'linewidth', 1.5, 'markersize',10)
% semilogy(s_irnt3(1:20)/s_irnt3(1), '+', 'linewidth', 1.5, 'markersize',10)
% semilogy(s_irnt4(1:20)/s_irnt4(1), 'd', 'linewidth', 1.5, 'markersize',10)
% legend('k = 1','k = 2','k = 3','k = 4')
% axis([1 20 10^-4.5 1])
% set(gca, 'fontsize',16), xticks([5 10 15 20]), yticks([10^-4,10^-2,1])


%% plot relative errors

figure
semilogy(0:length(info_gmres.Enrm),[1;info_gmres.Enrm],'linewidth',1.5)
hold on
semilogy(0:length(Enrm_tot_irn),[1;Enrm_tot_irn],'linewidth',1.5)
% semilogy(0:length(Enrm_tot_irnt),[1;Enrm_tot_irnt],'linewidth',1)
semilogy(0:length(Enrm_fnnrt),[1;Enrm_fnnrt],'linewidth',1.5)
% semilogy(0:length(Enrm_fnnrt),[1;Enrm_fnnr],'linewidth',1)
% semilogy(Enrm_fnnrt,'linewidth',1)
semilogy(0:length(Enrm_LRgm),[1;Enrm_LRgm],'linewidth',1.5)
semilogy(0:length(RelErr_RS),[1;RelErr_RS],'linewidth',1.5)
semilogy(0:length(RelErr_SVT),[1;RelErr_SVT],'linewidth',1.5)
% semilogy(RelErr_RS,'linewidth',1)

% [minE,ind] = min(info_gmres.Enrm);
% semilogy(ind,minE,'s')
% [minE,ind] = min(Enrm_tot_irn);
% semilogy(ind,minE,'s')
% [minE,ind] = min(Enrm_tot_irnt);
% semilogy(ind,minE,'s')
% [minE,ind] = min(info_gmres.Enrm);
% semilogy(ind,minE,'s')

legend('GMRES','IRN-GMRES-NNR','FGMRES-NNR','LR-FGMRES','RS-LR-GMRES','SVT')
axis([0 200 0.2 1])
xticks([0 50 100 150 200])
yticks([0.2 0.5 1])
set(gca,'fontsize',14)
set(gcf, 'Position',  [200, 200, 900, 300])
legend('location','southeastoutside')



%%

figure, imagesc(Xb_gm), axis image, axis off
% figure, mesh(Xb_gm(128-40:128+40, 128-40:128+40)), axis([0 80 0 80 -0.1 1])\

% Xb_irn = reshape(Xfinals_irn(:,1), n, n); 
figure, imagesc(Xb_lrfgm(128-25:128+25,128-25:128+25)), axis image, axis off
Xb_irn = reshape(Xfinals_irn(:,2), n, n); 
figure, imagesc(Xb_irn(128-25:128+25,128-25:128+25)), axis image, axis off
Xb_irn = reshape(Xfinals_irn(:,3), n, n); 
figure, imagesc(Xb_irn(128-25:128+25,128-25:128+25)), axis image, axis off
Xb_irn = reshape(Xfinals_irn(:,4), n, n); 
figure, imagesc(Xb_irn(128-25:128+25,128-25:128+25)), axis image, axis off

% figure, mesh(Xb_irn(128-40:128+40, 128-40:128+40)), axis([0 80 0 80 -0.1 1 ])
figure, imagesc(Xb_irn), axis image, axis off
% figure, mesh(Xb_irnt(128-40:128+40, 128-40:128+40)), axis([0 80 0 80 -0.1 1 ])
figure, imagesc(Xb_fnnr), axis image, axis off
% figure, mesh(Xb_fnnr(128-40:128+40, 128-40:128+40)), axis([0 80 0 80 -0.1 1 ])
figure, imagesc(Xb_lrfgm), axis image, axis off
figure, imagesc(Xb_rs), axis image, axis off

%%
% figure, mesh(Xb_gm(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
% set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
figure, mesh(Xb_irn(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
% figure, mesh(Xb_irnt(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
% set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
figure, mesh(Xb_fnnr(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
figure, mesh(Xb_lrfgm(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
figure, mesh(Xb_rs(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
figure, mesh(Xb_svt(128-band:128+band, 128-band:128+band)), axis([0 50 0 50 -0.1 1 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 1])
