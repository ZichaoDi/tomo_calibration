function [WRec, S_shiftback, drift, driftAll] = recTompographyWithDrift_rotation(S, paraTwist, WSz, driftGT)
%%TODO: comments
% TODO: merge recTompographyWithDrift_rotation() and recTompographyWithDrift()

% currently keep two sets of main calib/reconstruction algorithm and forward, backward model for beam drift and rotation error. 
% rotation error: recTompographyWithDrift_rotation(), shiftForward(), shiftBackward() 
% beam drift: recTompographyWithDrift(), interpSino(), calcDrift()
% TODO: merge thoese two sets, especially recTompographyWithDrift_rotation() and recTompographyWithDrift()

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
    S0_GT =  shiftForward(S, -driftGT);

end


maxIter = 10; nDisplay = 4;
%maxIter = 20; nDisplay = 4;
%maxIter = 1; nDisplay = 3;

L0 = XTM_Tensor_XH(WSz, NTheta, NTau); % foward model without drift
LNormalizer = full(max(sum(L0,2)));% Compute the normalizer 0.0278 for WSz=100*100, NTheta,NTau=45*152 for 1st scan of no drift so that maximum row sum is 1 rather than a too small number
L0 = L0/LNormalizer;


driftAll = zeros(NTheta, maxIter); 
WPrev = zeros(WSz);
figNo = 101; 
paraTwist2 = paraTwist;
drift = zeros(NTheta, 1);
upsampling_factor = 100; % Upsampling factor (integer). Images will be registered to   within 1/usfac of a pixel.
is_use_exact_forward_model = true; % true: 20.36dB,  false: 19.64dB for phantom 100x100, maxDrift =1

for k = 1:maxIter        % When use 10 Iter, wny after 10 iter, may slightly worse?
    %paraTwist2{6} = paraTwist{6}*1e2; %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
    
    % we use larger regularizer weight regWt*scaleBegin, and in the end (k=maxIter), we use smaller regularizer weight regWt*scaleEnd
    % This is the key for good reconstruction. 

    scaleBegin = 1e0; %1e2; %1e2 in the beginning(k=1), 
    scaleEnd = 1e0; %1e0~1e1 in the end (k=maxIter),
    idx = find(strcmp(paraTwist, 'regWt'))+1; % 6 etc.
    if maxIter > 1
        paraTwist2{idx} = paraTwist{idx}*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)); %{'xSz'}    {1×2 double}    {'regFun'}    {'TV'}    {'regWt'}    {[4.0000e-07]} 
    else
        paraTwist2{idx} = paraTwist{idx}*scaleEnd;
    end
       
    % (I) Update WRec, given drift    
    S_shiftback = shiftForward(S, -drift);    % shift back to recover undrifted sino        
    if is_use_exact_forward_model % or k > 5
        L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
        WRec = solveTwist(S(:), L, 'x0', WPrev, paraTwist2{:}); % 5 iters: 21.43/23.45dB without exact scan model, drift error 7.32dB        
    else        
        WRec = solveTwist(S_shiftback(:), L0, 'x0', WPrev, paraTwist2{:}); % 5 iters: 21.43/23.45dB without exact scan model, drift error 7.32dB        
    end
    
        
    % (II) Update drift, given WRec    
    S0_Rec = reshape(L0*WRec(:), SSz);
    drift =  shiftBackward(S0_Rec, S, upsampling_factor);
    driftAll(:, k) = drift;

	if rem(k,nDisplay)==1 || nDisplay==1 || k==maxIter
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, WRec}); multAxes(@title, {'Ground Truth', sprintf('Rec %dth, psnr=%.2fdB', k, difference(WRec,WGT))}); linkAxesXYZLimColorView();
        figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S0_GT, S_shiftback}); multAxes(@title, {'Sino-No-drift-GT', sprintf('Sina-Rec, psnr=%.2fdB', difference(S_shiftback,S0_GT))}); linkAxesXYZLimColorView(); 
        figure(figNo);figNo=figNo+1; plot(driftGT, 'b'); hold on; plot(drift, 'r'); legend('drift-GT', 'drift-Rec');
    end

    % compare current drift with that of previous iteration, and compare it with groundtruth
    if k>1
        idxNoZero = 1:NTheta; %idxNoZero = 25:125;
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end
    
    if (k>2) && norm(driftAll(:,k)-driftAll(:,k-1), 2) < 1/upsampling_factor
        break % converged to a drift, as our drift resolution is 
    end
    
    WPrev = WRec;
end


if k < size(driftAll, 2); driftAll = driftAll(:, 1:k); end % remove all the unused driftAll

% show final result
S_shiftback = shiftForward(S, -drift);    % shift back to recover undrifted sino    
figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, WRec}); multAxes(@title, {'Ground Truth', sprintf('Rec Final, psnr=%.2fdB', difference(WRec,WGT))}); linkAxesXYZLimColorView(); 
figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S0_GT, S_shiftback}); multAxes(@title, {'Sino-Measured', sprintf('Sina-Rec Final, psnr=%.2fdB', difference(S_shiftback,S0_GT))}); linkAxesXYZLimColorView(); 
figure(figNo);figNo=figNo+1; plot(driftGT, 'b'); hold on; plot(drift, 'r'); legend('drift-GT', 'drift-Rec');
tilefigs

