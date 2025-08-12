function [WRec, SRec, drift, driftAll,info] = recTompographyWithDrift_TypeIII(S, paraTwist, WSz, driftGT,maxIter,lambset,normeps,tolLamb,tolW)
%% TODO: comments
% TODO: merge recTompographyWithDrift_rotation() and recTompographyWithDrift()

% currently keep two sets of main calib/reconstruction algorithm and forward, backward model for beam drift and rotation error. 
% rotation error: recTompographyWithDrift_rotation(), shiftForward(), shiftBackward() 
% beam drift: recTompographyWithDrift(), interpSino(), calcDrift()
% TODO: merge thoese two sets, especially recTompographyWithDrift_rotation() and recTompographyWithDrift()
% 2020.05.15_fixed bug of reversed use of calcDriftTypeIII
% 2020.06.26_changed the way how many intermidiate result are shown as figures on screen
%            added total_show as the total number of intermediate iteration result shown on screen (not every iteration show)
%            and is_iters_show is the bool array of whether an iteration should be shown

[NTheta, NTau] = size(S); SSz = [NTheta, NTau];

% extract GT from paraTwist
idx = find(strcmp(paraTwist, 'xGT'));
if ~isempty(idx) && ~isempty(paraTwist{idx+1})
    WGT = paraTwist{idx+1};
    assert(isempty(WSz) || all(WSz == size(WGT)), 'specified WSz does not match the size of XGT in paraTwist');
else
    % WGT not defined, we must have WSz defined
    WGT = zeros(WSz);
end


if ~exist('driftGT', 'var') || isempty(driftGT)
    driftGT = zeros(NTheta, 1);
    S0_GT = zeros(size(S));
else        
    S0_GT = interpSinoTypeIII(S, -driftGT);
end

%maxIter = 10; nDisplay = 5;
%maxIter = 20; nDisplay = 4;
%maxIter = 1; nDisplay = 3;
% maxIter = 40; 
total_show = 2;
is_iters_show = false(maxIter, 1); 
is_iters_show(round(linspace(1, maxIter, total_show))) = true;


L0 = XTM_Tensor_XH(WSz, NTheta, NTau); % foward model without drift
LNormalizer = full(max(sum(L0,2)));% Compute the normalizer 0.0278 for WSz=100*100 and 50*50, NTheta,NTau=45*152 for 1st scan of no drift so that maximum row sum is 1 rather than a too small number
L0 = L0/LNormalizer;


driftAll = zeros(NTheta, NTau, maxIter); 
WPrev = zeros(WSz); WWPrev = WPrev;
figNo = 101; 
paraTwist2 = paraTwist;
drift = zeros(NTheta, NTau);

is_use_exact_forward_model = true; %% TODO: CHANGE-back to true

[Ny,Nx] = size(WGT);
psnr_history = zeros(maxIter,1);
W_history = zeros(Ny*Nx,maxIter);
S_history = zeros(NTheta*NTau,maxIter);
lambda_history = zeros(maxIter,1);
misfit_history = zeros(maxIter,1);

regWtbegin = paraTwist{6};

for k = 1:maxIter        % When use 10 Iter, wny after 10 iter, may slightly worse?
    %paraTwist2{6} = paraTwist{6}*1e2; %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
    
    % we use larger regularizer weight regWt*scaleBegin, and in the end (k=maxIter), we use smaller regularizer weight regWt*scaleEnd
    % This is the key for good reconstruction. 
    
    idx = find(strcmp(paraTwist, 'regWt'))+1; % 6 etc.
    if strcmp(lambset,'old')
        scaleBegin = 1e2; %1e2; %1e2 in the beginning(k=1), 
        scaleEnd = 1e-1; %1e0~1e1 in the end (k=maxIter),
        if maxIter > 1
            paraTwist2{idx} = (k-1)/(maxIter-1)*(scaleEnd - 1)*paraTwist{idx} + paraTwist{idx};
%             paraTwist2{idx} = paraTwist{idx}*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)); %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
        else
            paraTwist2{idx} = paraTwist{idx}*scaleEnd;
        end
    elseif strcmp(lambset,'new')
        if k == 1
            paraTwist2{6} = paraTwist{6};
        else
            paraTwist3 = paraTwist2;
            paraTwist3{6} = 1e-15;
            paraTwist3{16} = 1e-4;
            WW = solveTwist(S(:), L, 'x0', WWPrev, paraTwist3{:}); 
            alpha(k) = norm(S(:) - L*WW(:));
            beta(k) = norm(S(:) - L*WRec(:));   
            paraTwist2{6} = abs((normeps - alpha(k))/(beta(k) - alpha(k)))*paraTwist2{6};
            WWPrev = WW;
        end
    end
    regWt = paraTwist2{idx};
    
    if k > 1  
        tol = abs(regWt - lambda_history(k-1))/lambda_history(k-1)
        if tol < tolLamb
                W_history(:,k:end) = [];
                S_history(:,k:end) = [];
                psnr_history(k:end) = [];
                misfit_history(k:end) = [];
                lambda_history(k:end) = [];
            break
        end
    end
       
    % (I) Update WRec, given drift    
    if is_use_exact_forward_model % or k > 5
        L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
        WRec = solveTwist(S(:), L, 'x0', WPrev, paraTwist2{:}); % 5 iters: 21.43/23.45dB without exact scan model, drift error 7.32dB        
    else
        % won't work, as S is not uniform beamlet, not sampled on regular interval of tau axis.
        S_shiftback = interpSinoTypeIII(S, -drift);    % shift back to recover undrifted sino, only used for reconstruction method without exact forward model        
        WRec = solveTwist(S_shiftback(:), L0, 'x0', WPrev, paraTwist2{:}); % 5 iters: 21.43/23.45dB without exact scan model, drift error 7.32dB        
        warning('no L generated to show comparisoin of ')
    end
    
    
%     if k > 1  
%         tol = norm(WRec(:) - W_history(:,k-1))/norm(W_history(:,k-1))
%         if tol < tolW
%                 W_history(:,k:end) = [];
%                 S_history(:,k:end) = [];
%                 psnr_history(k:end) = [];
%                 misfit_history(k:end) = [];
%                 lambda_history(k:end) = [];
%             break
%         end
%     end    
    
    
    
        
    % (II) Update drift, given WRec. 
    % Compute the sinogram without drift from current WRec
    S0_Rec = reshape(L0*WRec(:), SSz);    
    
    method = 'new'; % 'new' with TV, 'old'
    if strcmp(method, 'old')
        % (II.1) Old, assume all drift are independent, no prior
        % compare with current sinagrom (no drift but from WRec not WGT) with measured sinagrom (has drift)    
        % drift = - calcDriftTypeIII(S0_Rec, S);  % was wrong using calcDriftTypeIII(S0_Rec, S)
        drift = calcDriftTypeIII(S, S0_Rec);  % NOTE: this is not as good as drift = - calcDriftTypeIII(S0_Rec, S); why?
    else 
        % (II.2) New, assume all drifts are like an image, use TV or sparsity prior
        %drift = - calcDriftTypeIII_TVonD(S0_Rec, S, -driftGT, regWt);  
        drift = calcDriftTypeIII_TVonD(S, S0_Rec, driftGT, regWt);
    end
    
    
    driftAll(:, :, k) = drift;
    
    if is_iters_show(k)  % old code:  if rem(k, nDisplay)==1 || nDisplay==1 || k==maxIter
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, WRec}); multAxes(@title, {'Ground Truth', sprintf('Rec %dth, psnr=%.2fdB', k, difference(WRec,WGT))}); linkAxesXYZLimColorView();

        SRec = reshape(L*WRec(:), NTheta, NTau); 
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRec}); multAxes(@title, {'Sinogram-Measured-w-Drift', sprintf('Sina-from rec, psnr=%.2fdB', difference(SRec,S))}); linkAxesXYZLimColorView(); 
        
        mask = (S~=0) & (abs(driftGT)<1); driftGT_masked = driftGT; driftGT_masked(~mask) = 0; drift_masked = drift; drift_masked(~mask) = 0;
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {driftGT_masked, drift_masked, drift_masked-driftGT_masked}); multAxes(@title, {'drift-GT', 'drift-Rec', sprintf('psnr=%.2fdB', difference(drift_masked, driftGT_masked))}); linkAxesXYZLimColorView(); multAxes(@colorbar);
    end

    % compare current drift with that of previous iteration, and compare it with groundtruth
    if k>1        
        mask = (S~=0) & (abs(driftGT)<1);
        driftCur = driftAll(:,:, k); 
        driftPrev = driftAll(:,:, k-1); 
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftPrev(mask), driftCur(mask)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftCur(mask), driftGT(mask)));
%         if norm(driftCur-driftPrev, 2) < NTheta*NTau*eps
%             break % converged to a drift, as our drift resolution is 
%         end
    end
    
    
    WPrev = WRec;
    W_history(:,k) = WRec(:);
    S_history(:,k) = L*WRec(:);
    psnr_history(k) = psnr(WRec,WGT);
    resid = S_history(:,k)-S(:);
    misfit_history(k) = 0.5*(resid'*resid);
    lambda_history(k) = regWt;
    
end

info.W_history = W_history;
info.S_history = S_history;
info.psnr_history = psnr_history;
info.misfit_history = misfit_history;
info.lambda_history = lambda_history;


if k < size(driftAll, 2); driftAll = driftAll(:, :, 1:k); end % remove all the unused driftAll
