function [W_history,S_history,info] = recTomoDrift_Chang(WGT,L0,S,maxIter,...
                            solver,maxDrift,LNormalizer,driftGT,params)
                        
% convert the script recTomographyWithDrift.m to a function

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end


[Ny,Nx] = size(WGT);
[NTheta, NTau] = size(S);
WSz = size(WGT);


driftAll = zeros(NTau, maxIter); 
L = L0; 

WPrev = zeros(Ny*Nx,1);
W0 = WPrev;

scaleBegin = 1e2; %1e2 in the beginning(k=1), 
scaleEnd = 1e0; %1e0~1e1 in the end (k=maxIter),

psnr_history = zeros(maxIter,1);


lambda_history = zeros(maxIter,1);
lambda_ind = 0;

W_history = zeros(Ny*Nx,maxIter);
S_history = zeros(NTheta*NTau,maxIter);
obj_history = zeros(maxIter,1);

P_history = cell(maxIter);

reg = params.reg;



for k = 1:maxIter 
    fprintf('iteration = %d ', k)

%     if maxIter > 1
%         params.reg = reg*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)) ;
%     else
%         params.reg = reg*scaleEnd;
%     end

    if strcmp(solver,'FLSQR-NNR')
        [WAll, EnrmAll,LambAll] = flsqr_nnr(L, S(:), WPrev, WGT(:), params);
    elseif strcmp(solver,'IRN-LSQR-NNR')
        [WAll, ~,~, EnrmAll,~,~,~,~,LambAll] = irn_lsqr_nnr(L, S(:), WGT(:), WPrev, params);
    elseif strcmp(solver,'LSQR')
        params.kappa = 0; params.kappaB = 0; params.iftrun = 0; params.tau = 0;
        [WAll, EnrmAll,~,~,LambAll] = LRlsqr(L, S(:), WPrev,WGT(:), params);
    end
    
    
    [~,ind] = min(EnrmAll);
    W = WAll(:,ind);    
    W = reshape(W,Ny,Nx);
    W(W<0) = 0;
    WPrev = W(:);

    
    obj_history(k) = norm(L*W(:)-S(:));
    
    lambda_history(k) = LambAll(k);
%     lambda_indnew = lambda_ind + length(LambAll);
%     lambda_history(lambda_ind + 1:lambda_indnew) = LambAll;
%     lambda_ind = lambda_indnew;
    
    fprintf('best solution given by iteration %d of inner solver\n', ind)
        
    S0 = reshape(L*W(:), NTheta, NTau);
    drift =  calcDrift(S, S0, maxDrift);
    driftAll(:, k) = drift;
    
    psnr_history(k) = psnr(W, WGT);
    W_history(:,k) = W(:);
    S_history(:,k) = L*W(:);

    if k < maxIter
        P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L; % update L^k
        P_history{k} = P;
        %L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
    end
    
    
    
    if k>1
        idxNoZero = 25:125;
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end 
    
%     if k > 1 && psnr_history(k) < psnr_history(k-1)
%         psnr_history(k:end) = [];
%         W_history(:,k:end) = [];
%         S_history(:,k:end) = [];
%         psnr_history(k:end) = [];
%         obj_history(k:end) = [];
%         lambda_history(k:end) = [];
%         break
%     end
    
end


% Final reconstruction use L1-norm
% P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L0;


% if strcmp(solver,'FLSQR-NNR')
%     [WAll,EnrmAll] = flsqr_nnr(L, S(:), W0, WGT(:), params);
% elseif strcmp(solver,'IRN-LSQR-NNR')
%     [WAll, ~,~, EnrmAll] = irn_lsqr_nnr(L, S(:), WGT(:), W0, params);
% end
% 
% [~,ind] = min(EnrmAll);
% W = reshape(WAll(:,ind),Ny,Nx);
% W(W<0) = 0;
% figure; imagesc(W); axis image; axis off
% title(sprintf('%s PSNR=%.2fdB', 'Use Interpolated Forward Model', difference(W, WGT)));
% W_history(:,end) = W(:);
% S_history(:,end) = L*W(:);

% lambda_history(lambda_ind+1:end) = [];

info.psnr_history = psnr_history;
info.lambda_history = lambda_history;
info.obj_history = obj_history;
info.P_history = P_history;

