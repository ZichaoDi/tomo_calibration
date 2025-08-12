function [WRec, SRec, L, drift, driftAll] = recTompographyWithDrift(S, paraTwist, WSz, maxDrift, driftGT)
% S:  Measured sinogram with drift
% paraTwist:  required
% WSz:  size of object, optional, but if paraTwist doesn't contain valid xGT, WSz must be defined
% maxDrift:   optional (default 1)
% driftGT: opitional (default all 0)

[NTheta, NTau] = size(S); SSz = [NTheta, NTau];

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end


% extract GT from paraTwist
idx = find(strcmp(paraTwist, 'xGT'));
if ~isempty(idx) && ~isempty(paraTwist{idx+1})
    WGT = paraTwist{idx+1};
end
isValid_WGT = exist('WGT', 'var') && ~isempty(WGT);
isValid_WSz = exist('WSz', 'var') && ~isempty(WSz);
% We must have either valid WGT (can infer WSz) or valid WSz (can set WGT to 0 if not exist)
if ~isValid_WGT && ~isValid_WSz    
    error('Both WGT and WSz not specified!')
elseif ~isValid_WGT && isValid_WSz
    WGT = zeros(WSz);
elseif isValid_WGT && ~isValid_WSz
    WSz = size(WGT);
else
    % both valid
    assert(all(WSz == size(WGT)), 'specified WSz does not match the size of XGT in paraTwist');
end

if ~exist('driftGT', 'var') || isempty(driftGT)
    driftGT = zeros(NTau, 1);
end


L0 = XTM_Tensor_XH(WSz, NTheta, NTau); % foward model without drift
LNormalizer = full(max(sum(L0,2)));% Compute the normalizer 0.0278 for WSz=100*100, NTheta,NTau=45*152 for 1st scan of no drift so that maximum row sum is 1 rather than a too small number
L0 = L0/LNormalizer;

maxIter = 10; nDisplay = 4;
%maxIter = 1; nDisplay = 3;

is_use_exact_forward_model = true;
driftAll = zeros(NTau, maxIter); 
WPrev = zeros(WSz);
figNo = 101; 
paraTwist2 = paraTwist;
L = L0;
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
    idx = find(strcmp(paraTwist, 'regWt'))+1; % 6 etc.
    if maxIter > 1
        paraTwist2{idx} = paraTwist{idx}*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)); %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
    else
        paraTwist2{idx} = paraTwist{idx}*scaleEnd;
    end
    
    % (I) Update WRec     
    %WRec = solveWendy(S, L, WGT, WPrev); %bad, should never use solveWendy    
    WRec = solveTwist(S(:), L, 'x0', WPrev, paraTwist2{:}); % 5 iters: 21.43/23.45dB without exact scan model, drift error 7.32dB    
    WPrev = WRec; 
% Test for different initialization
%     if k==1  
%         WGTMax = max(WGT(:)); WRec = WGTMax*imnoise(WGT/WGTMax, 'gaussian', 0, 0.01); %0.01*max(WGT(:))
%         WRec = WGT;
%     end
    if rem(k,nDisplay)==1 || nDisplay==1 || k==maxIter
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, WRec}); multAxes(@title, {'Ground Truth', sprintf('Rec %dth, psnr=%.2fdB', k, difference(WRec,WGT))}); linkAxesXYZLimColorView();
        SRec = reshape(L*WRec(:), SSz); 
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRec}); multAxes(@title, {'Sino-Measured', sprintf('Sina-Rec, psnr=%.2fdB', difference(SRec,S))}); linkAxesXYZLimColorView(); 
    end

    %%
    % (II) Update L
    % compute the drift vector as d+alpha, where d is the integer part, and 0 <= alpha <1 is the fractional part        
    S0 = reshape(L0*WRec(:), SSz);
    drift =  calcDrift(S, S0, maxDrift);
    driftAll(:, k) = drift;
    
    % update L^k given drift, compute the L matrix either from apporximation or exact
    if is_use_exact_forward_model % or k > 5
        L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
    else        
        P =  genDriftMatrix(drift, NTheta, NTau); 
        L =  P*L0; 
    end
    % compare current drift with that of previous iteration, and compare it with groundtruth
    if k>1 && all(size(drift) == size(driftGT)) % some debug use 2D driftGT        
        idxNoZero = 25:125;        
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end    

end


if k < size(driftAll, 2); driftAll = driftAll(:, 1:k); end % remove all the unused driftAll
    

%% show final result
if all(size(drift) == size(driftGT)) 
    figure(figNo);clf; figNo=figNo+1; plot([driftGT, drift]); legend('drift: GT', 'drift: Rec')
    SRec = reshape(L*WRec(:), SSz); %%SRec2 = interpSino(S0, drift, 'linear'); difference(SRec2, SRec) % 317.24dB matches well
    figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, WRec}); multAxes(@title, {'Ground Truth', sprintf('Rec Final, psnr=%.2fdB', difference(WRec,WGT))}); linkAxesXYZLimColorView(); 
    figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRec}); multAxes(@title, {'Sino-Measured', sprintf('Sina-Rec Final, psnr=%.2fdB', difference(SRec,S))}); linkAxesXYZLimColorView(); 
end
tilefigs

end
%%

%% Final reconstruction use L1-norm
function showResult(drift, NTheta, NTau, LAll, paraTwist, SAll, WGT)
P =  genDriftMatrix(drift, NTheta, NTau); L =  P*LAll{1};
regWt = paraTwist{find(strcmp(paraTwist, 'regWt'))+1};
regType = paraTwist{find(strcmp(paraTwist, 'regFun'))+1};
maxIterA = paraTwist{find(strcmp(paraTwist, 'maxIterA'))+1};
W2 = solveTwist(SAll(:,:,2), L, paraTwist{:}); figure(figNo+114); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Interpolated Forward Model', difference(W, WGT), regType, regWt, maxIterA));

% How good if we use exact model in final reconstruction, rather than use interpolation
LExactfromDrift = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
difference(LExactfromDrift*WGT(:), LAll{2}*WGT(:)); % psnr/psnrSparse = 42.41dB/35.70dB
W2 = solveTwist(SAll(:,:,2), LExactfromDrift, paraTwist{:}); figure(figNo+115); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Exact Forward Model from Calculated Drift', difference(W, WGT), regType, regWt, maxIterA));

end

