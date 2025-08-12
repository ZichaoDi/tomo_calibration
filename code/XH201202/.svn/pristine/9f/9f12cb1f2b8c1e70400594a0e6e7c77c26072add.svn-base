if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

L0 = LAll{1};
S = SAll(:,:,2); % observation, won't change
driftGT = driftGTAll{2};


[NTheta, NTau] = size(S);
M = NTheta*NTau; 
WSz = size(WGT);


%drift = (2*rand(NTau, 1)- 1) * 1;  P =  genDriftMatrix(drift, NTheta, NTau);  L =  P*L0;
maxIter = 20; nDisplay = 4;
%maxIter = 1; nDisplay = 3;
driftAll = zeros(NTau, maxIter); 
L = L0;  %P = sparse(1:M, 1:M, 1); % initialize as identity matrix

XPrev = zeros(Ny, Nx);
figNo = 101; 
paraTwist2 = paraTwist;

for k = 1:maxIter        % wny after 10 iter, may slightly worse?
    %paraTwist2{6} = paraTwist{6}*1e2; %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
    
    % we use larger regularizer weight regWt*scaleBegin, and in the end (k=maxIter), we use smaller regularizer weight regWt*scaleEnd
    % This is the key for good reconstruction. 
    % For comparison with maxIter=10, 
    % maxDrift = 1
    % scaleBegin/End = 1e2/1e0 gave 27.35dB but scaleBegin/End =1e0/1e0 only gave 23.16dB for final reconstruction of 100*100 brain
    % scaleBegin/End = 1e2/1e0 gave 21.26dB but scaleBegin/End =1e0/1e0 only gave 20.86dB for final reconstruction of 100*100 for phantom
    % scaleBegin/End = 1e2/1e1 gave 23.27dB but scaleBegin/End =1e1/1e1 !!!! gave 24.09dB for final reconstruction of 100*100 for phantom  (1e3/1e1 gave 20.09dB bad?1)
    % maxDrift = 3
    % scaleBegin/End = 1e2/1e0 gave 24.78dB but scaleBegin/End =1e0/1e0 only gave 14.36dB for final reconstruction of 100*100 brain
    
    scaleBegin = 1e2; %1e2 in the beginning(k=1), 
    scaleEnd = 1e0; %1e0~1e1 in the end (k=maxIter),
    if maxIter > 1
        paraTwist2{6} = paraTwist{6}*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)); %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
    else
        paraTwist2{6} = paraTwist{6}*scaleEnd;
    end
    % (I) Update X 
    if k <= 0
        X = solveWendy(S, L, WGT, XPrev); % 5 iters 21.12/23dB, 10 iters 21.26/23.45dB without exact scan model
    else
        X = solveTwist(S(:), L, 'x0', XPrev, paraTwist2{:}); % 5 iters: 21.43/23.45dB without exact scan model, drift error 7.32dB
    end
    XPrev = X;
    if k==1
        %WGTMax = max(WGT(:)); X = WGTMax*imnoise(WGT/WGTMax, 'gaussian', 0, 0.01); %0.01*max(WGT(:))
        %X = WGT;
    end  
    if rem(k,nDisplay)==1 || nDisplay==1
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, X}); multAxes(@title, {'Ground Truth', sprintf('Rec %dth, psnr=%.2fdB', k, difference(X,WGT))}); linkAxesXYZLimColorView();     
        SRec = reshape(L*X(:), NTheta, NTau); 
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRec}); multAxes(@title, {'Sinogram-Measured', sprintf('Sina-from rec, psnr=%.2fdB', difference(SRec,S))}); linkAxesXYZLimColorView(); 
    end

    %%
    % (II) Update L
    % compute the drift vector as d+alpha, where d is the integer part, and 0 <= alpha <1 is the fractional part        
    S0 = reshape(L0*X(:), NTheta, NTau);
    drift =  calcDrift(S, S0, maxDrift);
    driftAll(:, k) = drift;
    % given drift, compute the L matrix either from apporximation or exact
    if k <= 0%5
        P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L0; % update L^k
    else        
        L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
    end
    % compare driftOld and drift
    if k>1
        idxNoZero = 25:125;
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end    

end

%%
% remove all the unused driftAll
if k < size(driftAll, 2)
    driftAll(:, k) = [];
end


%% show final result
figure(figNo);clf; figNo=figNo+1; plot([driftGT, drift]); legend('drift: GT', 'drift: Rec')
%%
SRec = reshape(L*X(:), NTheta, NTau); %%SRec2 = interpSino(S0, drift, 'linear'); difference(SRec2, SRec) % 317.24dB matches well
figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, X}); multAxes(@title, {'Ground Truth', 'Reconstruction Final'}); linkAxesXYZLimColorView(); 
figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRec}); multAxes(@title, {'Sinogram-Measured', sprintf('Sina-from rec, psnr=%.2fdB', difference(SRec,S))}); linkAxesXYZLimColorView(); 
%figure, multAxes(@imagesc, {LAll{2}, L}); multAxes(@title, {'L Matrix - Ground Truth ', 'L Matrix - Computed'}); linkAxesXYZLimColorView();
%figure(figNo);figNo=figNo+1; multAxes(@imagesc, {LAll{2}, L, L0}); multAxes(@title, {'L Matrix - Ground Truth ', 'L Matrix - Computed w Drift', 'L0 - No Drift'}); linkAxesXYZLimColorView();
difference(S, SRec)
%difference(L,  LAll{2});


%% Final reconstruction use L1-norm
P =  genDriftMatrix(drift, NTheta, NTau); L =  P*LAll{1};
%difference(L, LAll{2})  % psnr/psnrSparse = 28.91dB/8.38dB
difference(L*WGT(:), LAll{2}*WGT(:)) % psnr/psnrSparse = 34.99dB/32.49dB
W2 = solveTwist(SAll(:,:,2), L, paraTwist{:}); figure(figNo+114); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Interpolated Forward Model', difference(W, WGT), regType, regWt, maxIterA));

% How good if we use exact model in final reconstruction, rather than use interpolation
LExactfromDrift = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
%difference(LExactfromDrift, LAll{2}); % psnr/psnrSparse = 32.53/12.00dB
difference(LExactfromDrift*WGT(:), LAll{2}*WGT(:)); % psnr/psnrSparse = 42.41dB/35.70dB
W2 = solveTwist(SAll(:,:,2), LExactfromDrift, paraTwist{:}); figure(figNo+115); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Exact Forward Model from Calculated Drift', difference(W, WGT), regType, regWt, maxIterA));
%W2 = solveWendy(SAll(:,:,2), LExactfromDrift, WGT); figure(figNo+15); W = W2;  imagesc(W); axis image; title(sprintf('Rec Wendy %s, PSNR=%.2fdB', 'Use Exact Forward Model from Drift', difference(W, WGT)));

%%
paraTwist2 = paraTwist;
paraTwist2{6} = paraTwist2{6}*10;
W2 = solveTwist(SAll(:,:,2), LExactfromDrift, paraTwist2{:}); figure(figNo+116); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Exact Forward Model from Calculated Drift', difference(W, WGT), regType, paraTwist2{6}, maxIterA));


%% GT
%PGT =  genDriftMatrix(driftGT, NTheta, NTau);        
%LGT =  PGT*L0; % update L^k
%SRecGT = reshape(LGT*WGT(:), NTheta, NTau); 
%figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRecGT}); multAxes(@title, {'Sinogram-Measured', sprintf('Sina-from recGT, psnr=%.2fdB', difference(SRecGT,S))}); linkAxesXYZLimColorView(); 
%figure, multAxes(@imagesc, {LAll{2}, L}); linkAxesXYZLimColorView;
%figure, multAxes(@imagesc, {reshape(LAll{2}*X(:), NTheta, []), reshape(Lrec*X(:), NTheta, [])}); linkAxesXYZLimColorView;
% W = solveWendy(SAll(:,:,2), Lrec, WGT);
% figure(222); imagesc(W); axis image; 

tilefigs

%%

