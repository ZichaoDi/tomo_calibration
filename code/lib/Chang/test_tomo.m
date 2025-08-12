 clear, clc
 rng(52);

% 
optomo = PRtomo('defaults');
optomo.phantomImage = 'smooth';
% optomo.angles = 0:2:179;
optomo.angles = 0:1:90;

n = 256; N = n^2;

[A, b, x, ProbInfo] = PRtomo(n, optomo);
% [~, ~, x] = PRtomo(n, optomo);
% 
% % LR solution
% Xex = reshape(x, n, n);
% [UX, SX, VX] = svd(Xex);
% tx = 50;
% Xext = UX(:,1:tx)*(SX(1:tx,1:tx)*VX(:,1:tx)');
% 
% optomo.phantomImage = Xext;
% [A, b, x, ProbInfo] = PRtomo(n, optomo);

%
nl = 1e-2;
bn = PRnoise(b,nl);

% p = 0.75;

reg = 0;
p2 = 0.8;

%%

opth.RegParam = reg;
opth.NoiseLevel = nl;
% opth.MaxIter = 100;
opth.NoStop = 'on';
opth.x_true = x;
opth.DecompOut = 'on';

% standard LSQR
nit = 100;
[X_lsqr, info_lsqr] = IRhybrid_lsqr(A, bn, 1:nit, opth);
% checking the rank 
rkX_lsqr = zeros(nit,1);
for i = 1:nit
    Xtemp = reshape(X_lsqr(:,i), n, n);
    sX = svd(Xtemp);
    rkX_lsqr(i) = length(sX(sX>1e-10));
end

%%
rktrunc = 10;
k = 100; eta=1.01; x0 = zeros(N,1);
parameters = struct('maxIt',k','kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
[X_flsqr, Enrm_flsqr, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% parameters = struct('maxIt',k','kappaB',rktrunc,'kappa',n,'iftrun',1,'tau',0,'reg',0,'eta',eta,'nl',nl);
% [X_flsqrvar, Enrm_flsqrvar, Rnrm_flsqrvar, rkX_flsqrvar] = LRlsqr(A, bn, x0, x, parameters);

%%
% rktrunc = 10;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k','kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr2, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);

% rktrunc = 8;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k','kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr3, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);
% 
% rktrunc = 10;
% k = 100; eta=1.01; x0 = zeros(N,1);
% parameters = struct('maxIt',k','kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg',reg,'eta',eta,'nl',nl);
% [X_flsqr, Enrm_flsqr4, Rnrm_flsqr, rkX_flsqr] = LRlsqr(A, bn, x0, x, parameters);



%%
% figure, semilogy(Enrm_flsqr)
% hold on
% semilogy(Enrm_flsqr2)
% semilogy(Enrm_flsqr3)
% semilogy(Enrm_flsqr4)
% legend('1','2','3','4')


%%
optnnr.p = 1;
optnnr.maxIt = 100;
optnnr.reg = reg;
optnnr.eta = eta;
optnnr.nl = nl;
optnnr.svdbasis = 'v';
optnnr.thr = 1e-3;
optnnr.gamma = 10^-10;
[X_flsqrnnr1, Enrm_flsqrnnr1, Rnrm_flsqrnnr1] = flsqr_nnr(A, bn, x0, x, optnnr);
%%
optnnr.p = p2;
[X_flsqrnnr1p, Enrm_flsqrnnr1p, Rnrm_flsqrnnr1p] = flsqr_nnr(A, bn, x0, x, optnnr);
%%
optnnr.p = 1;
optnnr.svdbasis = 'x';
[X_flsqrnnr5, Enrm_flsqrnnr5, Rnrm_flsqrnnr5] = flsqr_nnr(A, bn, x0, x, optnnr);
optnnr.p = p2;
[X_flsqrnnr5p, Enrm_flsqrnnr5p, Rnrm_flsqrnnr5p] = flsqr_nnr(A, bn, x0, x, optnnr);
% % optnnr.svdbasis = 6;
% % [X_flsqrnnr6, Enrm_flsqrnnr6, Rnrm_flsqrnnr6] = flsqr_nnr(A, bn, x0, x, optnnr);
% 
%

%%
optnnr.gamma = 10^-10;
optnnr.maxIt = 25;
optnnr.cycles = 4;
optnnr.thr = 1e-3;
optnnr.p = 1;
optnnr.reg = 0;
optnnr.weigthtype = 'sqrt'; 
optnnr.thrstop = 1e-8;
[X_irnlsqr_tot, X_irnlsqr_finals, X_irnlsqr_best, Enrm_irnlsqr, Rnrm_irnlsqr, ...
    cyclesIt, sXcycles, cyclesBest, Lambda_tot] = irn_lsqr_nnr(A, bn, x, x0, optnnr);
optnnr.p = p2;
[X_irnlsqrp_tot, X_irnlsqrp_finals, X_irnlsqrp_best, Enrm_irnlsqrp, Rnrm_irnlsqrp, ...
    cyclesItp, sXcyclesp, cyclesBestp, Lambda_totp] = irn_lsqr_nnr(A, bn, x, x0, optnnr);


%% other methods
[X_SVT,RelErr_SVT,RelRes_SVT] = SVT(A,bn,x,1,100,100,0.00008);

%% RS-LR-GMRES on normal equations
A_ne = @(x) A'*(A*x);
bn_ne = A'*bn;
tau = 10;
[X_RS, RelRes_RS, RelErr_RS] = RS_GMRES_LRP(A_ne,bn_ne,10,10,x0,tau,tau,x);
figure,semilogy(RelErr_RS)
min(RelErr_RS)
%%
%%%% COMSIDERING THE DISCREPANCY PRINCIPLE
% opth.RegParam = 'discrep';
% [X_lsqrq, info_lsqrq] = IRhybrid_lsqr(A, bn, 1:nit, opth);
% parameters = struct('maxIt',k','kappaB',rktrunc,'kappa',rktrunc,'iftrun',1,'tau',0,'reg','discrep','eta',eta,'nl',nl);
% [X_flsqrq, Enrm_flsqrq, Rnrm_flsqrq, rkX_flsqrq] = LRlsqr(A, bn, x0, x, parameters);
% optnnr.maxIt = 100;
% optnnr.svdbasis = 't';
% optnnr.p = 1;
% optnnr.reg = 'discrep';
% [X_flsqrnnr1q, Enrm_flsqrnnr1q, Rnrm_flsqrnnr1q] = flsqr_nnr(A, bn, x0, x, optnnr);
% optnnr.p = 0.75;
% [X_flsqrnnr1pq, Enrm_flsqrnnr1pq, Rnrm_flsqrnnr1pq] = flsqr_nnr(A, bn, x0, x, optnnr);
% optnnr.svdbasis = 'x';
% optnnr.p = 1;
% [X_flsqrnnr5q, Enrm_flsqrnnr5q, Rnrm_flsqrnnr5q] = flsqr_nnr(A, bn, x0, x, optnnr);
% optnnr.p=0.75;
% [X_flsqrnnr5pq, Enrm_flsqrnnr5pq, Rnrm_flsqrnnr5pq] = flsqr_nnr(A, bn, x0, x, optnnr);
% optnnr.maxIt = 25;
% optnnr.p = 1;
% [X_irnlsqrq_tot,X_irnlsqrq_finals, X_irnlsqrq_best, Enrm_irnlsqrq, Rnrm_irnlsqrq, ...
%     cyclesItq, sXcyclesq, cyclesBestq, Lambda_totq] = irn_lsqr_nnr(A, bn, x, x0, optnnr);
% optnnr.p = 0.75;
% [X_irnlsqrpq_tot,X_irnlsqrpq_finals, X_irnlsqrpq_best, Enrm_irnlsqrpq, Rnrm_irnlsqrpq, ...
%     cyclesItpq, sXcyclespq, cyclesBestpq, Lambda_totpq] = irn_lsqr_nnr(A, bn, x, x0, optnnr);
% 
%
% % % new version of NNR-FLSQR
% % parameters = struct('maxIt',k','svdbasis',1,'reg',0,'eta',eta,'nl',nl,'p',1);
% % tic;[X1, relerr1] = flsqr_nnr(A, bn, x0,x,parameters);toc;
% % parameters = struct('maxIt',k','svdbasis',2,'reg',0,'eta',eta,'nl',nl,'p',1);
% % tic;[X2, relerr2] = flsqr_nnr(A, bn, x0,x,parameters);toc;
% % parameters = struct('maxIt',k','svdbasis',3,'reg',0,'eta',eta,'nl',nl,'p',1);
% % tic;[X3, relerr3] = flsqr_nnr(A, bn, x0,x,parameters);toc;
% % parameters = struct('maxIt',k','svdbasis',4,'reg',0,'eta',eta,'nl',nl,'p',1);
% % tic;[X4, relerr4] = flsqr_nnr(A, bn, x0,x,parameters);toc;
% % parameters = struct('maxIt',k','svdbasis',5,'reg',0,'eta',eta,'nl',nl,'p',1);
% % tic;[X5, relerr5] = flsqr_nnr(A, bn, x0,x,parameters);toc;
% % parameters = struct('maxIt',k','svdbasis',6,'reg',0,'eta',eta,'nl',nl,'p',1);
% % tic;[X6, relerr6] = flsqr_nnr(A, bn, x0,x,parameters);toc;
% 

%%
%%%% plotting the results
X = reshape(x, n, n);
XX = insertShape(X,'rectangle',[1 1 64 64],'linewidth',1,'opacity',0.3);
figure, imagesc(rgb2gray(XX)), axis image, axis off
%figure, mesh(X), axis image, axis off


Xlsqr = reshape(info_lsqr.BestReg.X, n, n);
figure, imagesc(Xlsqr), axis image, axis off
%figure, mesh(Xlsqr), axis image, axis off


[~, ind] = min(Enrm_flsqr);
Xflsqr = reshape(X_flsqr(:,ind), n, n);
figure, imagesc(Xflsqr), axis image, axis off
%figure, mesh(Xflsqr), axis image, axis off

[~, ind] = min(Enrm_flsqrnnr1);
Xfnnr = reshape(X_flsqrnnr1(:,ind),n,n);
figure, imagesc(Xfnnr), axis image, axis off
%figure, mesh(Xfnnr), axis image, axis off

[~, ind] = min(Enrm_flsqrnnr1p);
Xfnnrp = reshape(X_flsqrnnr1p(:,ind),n,n);
figure, imagesc(Xfnnrp), axis image, axis off
%figure, mesh(Xfnnrp), axis image, axis off

[~, ind] = min(Enrm_irnlsqr);
Xirnlsqr = reshape(X_irnlsqr_tot(:,ind),n,n);
figure, imagesc(Xirnlsqr), axis image, axis off
%figure, mesh(Xirnlsqr), axis image, axis off

[~, ind] = min(RelErr_RS);
Xrs = reshape(X_RS(:,ind),n,n);
figure, imagesc(Xrs), axis image, axis off


% figure, mesh(X(1:64,1:64)), axis image, axis off
% figure, mesh(Xlsqr(1:64,1:64)), axis image, axis off
% figure, mesh(Xflsqr(1:64,1:64)), axis image, axis off
% figure, mesh(Xfnnr(1:64,1:64)), axis image, axis off
% figure, mesh(Xirnlsqr(1:64,1:64)), axis image, axis off

%%
cut = 64;
figure, mesh(X(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])
figure, mesh(Xlsqr(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])
figure, mesh(Xflsqr(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])
figure, mesh(Xfnnr(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])
figure, mesh(Xfnnrp(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])
figure, mesh(Xirnlsqr(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])
figure, mesh(Xrs(1:cut,1:cut)), axis([0 70 0 70 -0.1 0.35 ])
set(gca,'XTick',[],'YTick',[],'fontsize',14), zticks([-0.1 0.35])

%%
figure, 
semilogy(0:length(Enrm_flsqrnnr5),[1;Enrm_flsqrnnr5],'linewidth',1.5)
hold on

semilogy(0:length(info_lsqr.Enrm),[1;info_lsqr.Enrm],'linewidth',1.5)
% semilogy(Enrm_flsqrnnr5)
semilogy(0:length(Enrm_flsqr),[1;Enrm_flsqr],'linewidth',1.5)
semilogy(0:length(RelErr_RS),[1;RelErr_RS],'linewidth',1.5)
semilogy(0:length(RelErr_SVT),[1;RelErr_SVT],'linewidth',1.5)
legend('FLSQR-NNR','LSQR', 'LR-FLSQR' ,'RS-LR-GMRES','SVT')
legend('location','southeastoutside')
set(gca,'fontsize',14)
set(gcf, 'Position',  [200, 200, 900, 300])
axis([0 100 0.05 1])
yticks([0.05 0.1 0.5 1])

figure, 
semilogy(0:length(Enrm_flsqrnnr5),[1;Enrm_flsqrnnr5],'linewidth',1.5)

%semilogy(0:length(info_lsqr.Enrm),[1;info_lsqr.Enrm],'linewidth',1.5)
hold on
% semilogy(0:length(Enrm_flsqr),[1;Enrm_flsqr],'linewidth',1.5)
% semilogy(0:length(Enrm_flsqrnnr5),[1;Enrm_flsqrnnr5],'linewidth',1.5)
semilogy(0:length(Enrm_flsqrnnr1),[1;Enrm_flsqrnnr1],'linewidth',1.5)
semilogy(0:length(Enrm_flsqrnnr1p),[1;Enrm_flsqrnnr1p],'linewidth',1.5)
% semilogy(Enrm_flsqrnnr5,'linewidth',1)
% semilogy(Enrm_flsqrnnr5p,'linewidth',1)
semilogy(0:length(Enrm_irnlsqr),[1;Enrm_irnlsqr],'linewidth',1.5)
semilogy(0:length(Enrm_irnlsqrp),[1;Enrm_irnlsqrp],'linewidth',1.5)
%semilogy(Enrm_irnlsqr_t(1:100))
legend('FLSQR-NNR', 'FLSQR-NNR(v)', 'FLSQR-NNRp(v)', 'IRN-LSQR-NNR' ,'IRN-LSQR-NNRp')
legend('location','southeastoutside')
set(gca,'fontsize',14)
set(gcf, 'Position',  [200, 200, 900, 300])
axis([0 100 0.05 1])
yticks([0.05 0.1 0.5 1])
