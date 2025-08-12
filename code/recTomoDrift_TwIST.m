function [W_history,S_history,info] = recTomoDrift_TwIST(WGT,L0,L1,S,maxIter,solver,...
                            maxDrift,lambSet,regWtBegin,regWtEnd,driftGT,params,optnnr)
                        
% convert the script recTomographyWithDrift.m to a function

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end



[Ny,Nx] = size(WGT);
[NTheta, NTau] = size(S);
% WSz = size(WGT);

driftAll = zeros(NTau, maxIter); 
L = L0; 

W0 = zeros(Ny, Nx); WPrev = W0; WWPrev = W0;
paraTwist2 = params;

psnr_history = zeros(maxIter,1);
W_history = zeros(Ny*Nx,maxIter);
S_history = zeros(NTheta*NTau,maxIter);
lambda_history = zeros(maxIter,1);

regFun_history = zeros(maxIter,1);
misfit_history = zeros(maxIter,1);
etaeps_history = zeros(maxIter,1);

alpha = zeros(maxIter,1);
beta = zeros(maxIter,1);
lambda_final = zeros(maxIter+1,1);


regFun = @(f) TVnorm3D(f, [], 'isotropic');
misfitFun = @(resid) 0.5*(resid(:)'*resid(:));


if strcmp(solver,'FLSQR') == 1
    [WAll, EnrmAll,RnrmAll,~,lambdahis,alphahis] = flsqr_nnr(L, S(:), WPrev(:), WGT(:), optnnr);
    [~,ind] = min(EnrmAll);
    W = WAll(:,end);    
    W = reshape(W,Ny,Nx);
    W(W<0) = 0;
    alpha(1) = alphahis(ind);
    beta(1) = RnrmAll(ind);
    final_lambda(1) = lambdahis(end);
elseif strcmp(solver,'TwIST') == 1
	W = solveTwist(S(:), L, 'x0', WPrev, paraTwist2{:}); 
else
    [W,info] = IRhtv(L,S(:),1:50,struct('x0',WPrev(:),'x_true',WGT(:),'MaxIter',50,...
        'MaxIterIn',30,'MaxIterOut',4,'RegParam0',regWtBegin,'adaptConstr','tvnn',...
        'NoiseLevel',norm(L1*WGT(:)-S(:))/norm(S(:)),'RegParam','discrep','inSolver','lsqr'));
    [~,ind]=min(info.Enrm);
    W = W(:,ind);
    W = reshape(W,Ny,Nx);
    W(W<0) = 0;
end

W00 = W;

psnr0 = psnr(W,WGT);
lambda0 = paraTwist2{6};
regFun0 = regFun(W);
misFit0 = misfitFun(L*W(:)-S(:));

k = 1;
while 1
    
    S0 = reshape(L*W(:), NTheta, NTau);
    drift =  calcDrift(S, S0, maxDrift);
    driftAll(:, k) = drift;
    
    

    P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L; % update L^k
%     L = XTM_Tensor_XH(WSz, NTheta, NTau, drift);
%     LNormalizer = full(max(sum(L,2)));
%     L = L/LNormalizer;
    
    
    
    etaeps_history(k) = norm(L1*WGT(:)-S(:));    

    if strcmp(lambSet,'new') == 1
        if strcmp(solver,'TwIST') == 1
            paraTwist3 = paraTwist2;
            paraTwist3{6} = 1e-15;
            paraTwist3{16} = 1e-4;
            WW = solveTwist(S(:), L, 'x0', WWPrev, paraTwist3{:}); 
            alpha(k) = norm(S(:) - L*WW(:));        
            beta(k) = norm(S(:) - L*W(:));        
            paraTwist2{6} = abs((etaeps_history(k) - alpha(k))/(beta(k) - alpha(k)))*paraTwist2{6};
%             if paraTwistnew >= 0.05*paraTwist2{6}
%                 paraTwist2{6} = paraTwistnew;
%             end
            WWPrev = WW;
        else
            paraTwist2{6} = abs((etaeps_history(k)/norm(S(:)) - alpha(k))/(beta(k) - alpha(k)))*paraTwist2{6};
            optnnr.reg = paraTwist2{6};
        end
    end   
    
    
    if strcmp(lambSet,'old') == 1
        paraTwist2{6} = (k-1)/(maxIter-1)*(regWtEnd - regWtBegin) + regWtBegin;
        optnnr.reg = paraTwist2{6};
    elseif strcmp(lambSet,'fixed') == 1 
        paraTwist2{6} = regWtEnd;
    end
    lambda_history(k) = paraTwist2{6};

    
    WPrev = W;

   

    if strcmp(solver,'FLSQR') == 1
        if strcmp(lambSet,'discrep') == 1
            optnnr.reg = 'discrep';
        end
                
        [WAll, EnrmAll,RnrmAll,objhis,lambdahis,alphahis] = flsqr_nnr(L, S(:), WPrev(:), WGT(:), optnnr);
        [~,ind] = min(EnrmAll);
        W = WAll(:,end);    
        W = reshape(W,Ny,Nx);
        W(W<0) = 0;

        figure(k),plot(EnrmAll),title('relerr')            
        lambda_final(k+1) = lambdahis(end);
        
    

        if k < maxIter
            alpha(k+1) = alphahis(ind);
            beta(k+1) = RnrmAll(ind);
        end
    elseif strcmp(solver,'TwIST') == 1
        W = solveTwist(S(:), L, 'x0', WPrev, paraTwist2{:}); 
    else
        if strcmp(lambSet,'discrep') == 1
            [W,info] = IRhybrid_lsqr(L,S(:),1:50,struct('x0',WPrev(:),'x_true',WGT(:),'MaxIter',50,...
                'MaxIterIn',30,'MaxIterOut',4,'RegPram0',regWtBegin,'adaptConstr','tvnn',...
                'NoiseLevel', etaeps_history(k)/norm(S(:)),'RegParam','discrep','inSolver','lsqr','NoStop','on'));
            info.its
        else
            [W,info] = IRhtv(L,S(:),1:50,struct('x0',WPrev(:),'x_true',WGT(:),'MaxIter',50,...
                'MaxIterIn',30,'MaxIterOut',4,'RegPram0',regWtBegin,'adaptConstr','tvnn',...
                'NoiseLevel', etaeps_history(k)/norm(S(:)),'RegParam',paraTwist2{6},'inSolver','lsqr'));
        end

        [~,ind] = min(info.Enrm)
        W = W(:,ind);
    %     W(W<0) = 0;
        W = reshape(W,Ny,Nx);

        figure(k),plot(info.Enrm)
        figure(k+20),plot(info.Rnrm)
        figure(k+40),plot(info.RegP)
    
    end

    regFun_history(k) = regFun(W);
    resid = S(:)-L*W(:); %nres = norm(resid);
    misfit_history(k) = misfitFun(resid);
    
%     alpha = misfit_history(k) + lambda_history(k)*regFun_history(k);
%     beta = regFun_history(k);
%     paraTwist2{6} = abs((etaeps_history(k) - alpha)/beta)/k^(1/2);

   
    
    psnr_history(k) = psnr(W, WGT);
    W_history(:,k) = W(:);
    S_history(:,k) = L*W(:);
            

    if k == 1
        diffW(1) = norm(W(:)-W00(:));
    else
        diffW(k) = norm(W(:)-WPrev(:));
    end
    
    
    if k>1
        idxNoZero = 25:125;
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end   
%    
%     if k > 1 && psnr_history(k) < psnr_history(k-1)
%         psnr_history(k:end) = [];
%         W_history(:,k:end) = [];
%         S_history(:,k:end) = [];
%         break
%     end

    if k == 1
        criterion = abs(misfit_history(1) - misFit0)/misFit0
    else
        criterion = abs(misfit_history(k) - misfit_history(k-1))/misfit_history(k-1)
    end
    
    if k == maxIter
%     if criterion < 1e-3 || k == maxIter
        break
    end
    
    

    
    k = k + 1;
end

if k < maxIter
    psnr_history(k+1:end) = [];
    lambda_history(k+1:end) = [];
    regFun_history(k+1:end) = [];
    misfit_history(k+1:end) = [];
    etaeps_history(k+1:end) = [];
    alpha(k+1:end) = [];
    beta(k+1:end) = [];
    
    W_history(:,k+1:end) = [];
    S_history(:,k+1:end) = [];
end

info.psnr_history = [psnr0;psnr_history];
info.lambda_history = [lambda0;lambda_history];
info.regFun_history = [regFun0;regFun_history];
info.misfit_history = [misFit0;misfit_history];

info.etaeps_history = etaeps_history;
info.alpha = alpha;
info.beta = beta;
info.lambda_final = lambda_final;
info.driftAll = driftAll;
