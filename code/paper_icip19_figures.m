
%% load results and generate figures and tables, refer to paper_icip19.m which generate and save results
isShowFigs = 0; % 1/0
isUsedBestWtUncalib = 0; % 1/0
%retFile = 'ret_20190122_0056_NTheta45_3SmallNoiseOnSino_4Drift.mat';
retFile = 'ret_20190124_2020_NTheta45_4Drift_3SmallNoiseOnSino.mat'; % similar to ret_20190122_0056_NTheta45_3SmallNoiseOnSino
% retFile = 'ret_20190124_2339_NTheta45_4Drift_3SmallNoiseOnSino.mat'; % similar to ret_20190122_0056_NTheta45_3SmallNoiseOnSino
% retFile = 'ret_20190125_0236_NTheta45_4Drift_3SmallNoiseOnSino_BaselineBestregWt_Ours1_0e-5.mat'; % 
% retFile = 'ret_20190125_0227_NTheta45_4Drift_3SmallNoiseOnSino_BaselineBestregWt_Ours1_3e-5.mat'; % use this one 1_3e-5 not 1.0e-5
% retFile = 'ret_20190125_1305.mat';
% retFile = 'ret_20190125_1320.mat';

opengl('save','software'); %hack fix for crashed matlab due to the graphics driver; https://www.mathworks.com/matlabcentral/answers/157894-resolving-low-level-graphics-issues-in-matlab
% We compare performance in both peak signal to noise ratio (PSNR) and structural similarity (SSIM).
% In Figure 1, we shows the result of noiseless case, we  various amount of drift error: $d_{max} = 0.5\Delta, \Delta, 2\Delta, 3\Delta, 5\Delta$. 

% save(['ret_' curTime '.mat'], ...
%     'sampleNameAll', 'gaussianSTDAll', 'maxDriftAll', 'noise', ...
%     'WGTAll', 'SNoDriftAll', 'SDriftAll', 'driftGTAll', 'driftRecAll',  'WUncalibAll', 'WRecAll', 'SRecAll', ...
%     'maxIter', 'scaleBegin', 'scaleEnd'), 
load(retFile);

if ~exist('gaussianSTDAll', 'var') || isempty(gaussianSTDAll)
    gaussianSTDAll = [0];
end


NSample = numel(sampleNameAll);
NNoise = numel(gaussianSTDAll);
NMaxDrift = numel(maxDriftAll);


%% (1) Calculate Quantatitive PSNR and SSIM
% for each sample, plot the PSNR and drift error
mat0 = zeros(NSample, NMaxDrift, NNoise);
PSNR_RecAll = mat0;
SSIM_RecAll = mat0;
PSNR_UncalibAll = mat0;
SSIM_UncalibAll = mat0;
for nSample = 1:NSample
    sampleName = sampleNameAll{nSample};
    WGT = WGTAll{nSample};
    for nMaxDrift = 1:NMaxDrift
        maxDrift = maxDriftAll(nMaxDrift);                
        for nNoise = 1:NNoise
            gaussianSTD = gaussianSTDAll(nNoise);        
            % Reconstruction
            % drift = driftRecAll{nSample, nMaxDrift, nNoise};
            if isUsedBestWtUncalib
                WUncalib = WUncalibBestWtAll{nSample, nMaxDrift, nNoise}; % best weight
            else
                WUncalib = WUncalibAll{nSample, nMaxDrift, nNoise}; % fixed weight
            end
            WRec = WRecAll{nSample, nMaxDrift, nNoise} ;
            % SRec = SRecAll{nSample, nMaxDrift, nNoise};
            
            if isShowFigs
                figure; multAxes(@imagesc,{WRec, WUncalib(:,:,1), WUncalib(:,:,2), WUncalib(:,:,3)});
            end
            
            PSNR_RecAll(nSample, nMaxDrift, nNoise) = psnr(WRec, WGT);  % here we assume max(WGT(:)) == 1
            SSIM_RecAll(nSample, nMaxDrift, nNoise) = ssim(WRec, WGT);
            PSNR_UncalibAll(nSample, nMaxDrift, nNoise) =  psnr(WUncalib(:,:,3), WGT);
            SSIM_UncalibAll(nSample, nMaxDrift, nNoise) =  ssim(WUncalib(:,:,3), WGT);
        end
    end
end
% squeeze especially for NNoise == 1;
PSNR_RecAll = squeeze(PSNR_RecAll);
SSIM_RecAll = squeeze(SSIM_RecAll);
PSNR_UncalibAll = squeeze(PSNR_UncalibAll);
SSIM_UncalibAll = squeeze(SSIM_UncalibAll);
if isShowFigs
    tilefigs([], [], NSample*NNoise, NMaxDrift)
end

%% Print to format that suitable for table
tableCell = {};
for n = 1:2
    tableCell = [tableCell, {PSNR_RecAll(n,:,:), SSIM_RecAll(n,:,:), PSNR_UncalibAll(n,:,:), SSIM_UncalibAll(n,:,:)}];
end

%% each line of table is 12.34 for PSNR or .1234 for SSIM
tableMatrix = [];
for n = 1:numel(tableCell)
    curLine = squeeze(tableCell{n})';
    curLine = curLine(:)';
    tableMatrix = [tableMatrix; curLine];
    sep = ' & ';
    for m = 1:numel(curLine)        
        if mod(n,2) == 1
            fprintf('%.2f', curLine(m))
        else
            fprintf('%.4f', curLine(m))
        end
        if m < numel(curLine)
            fprintf(sep);
        end
    end
    fprintf('\\\\ \n');
end
PSNRDiff = tableMatrix([1 5], :) - tableMatrix([3 7], :);
SSIMDiff = tableMatrix([2 6], :) - tableMatrix([4 8], :);
mean(PSNRDiff(:))
mean(SSIMDiff(:))
tableMatrix

%% Save figure, we have NSample*NNoise*NMaxDrift result, only save selected noise level and maxDrift lever
retFolder = fullfile('..', '..', 'writing', 'driftedTomography', 'figs'); %'xx\Dropbox\Tomography\writing\driftedTomography\figs';
driftSelect = [1, NMaxDrift]; % NMaxDrift = 4; use 1:NMaxDrift to select all.
noiseSelect = [1:NNoise];  % [1, NNoise] NNoise = 3, use 1:NNoise to select all.

% only work on selected subset of noise and drift levels
maxDriftSelect = maxDriftAll(driftSelect);
gaussianSTDAllSelect = gaussianSTDAll(noiseSelect);
WUncalibAllSelect = WUncalibAll(:, driftSelect, noiseSelect);
if isUsedBestWtUncalib
    WUncalibBestWtAllSelect = WUncalibBestWtAll(:, driftSelect, noiseSelect);
end

WRecAllSelect = WRecAll(:, driftSelect, noiseSelect);
SDriftAllSelect = SDriftAll(:, driftSelect, noiseSelect);
SNoDriftAllSelect = SNoDriftAll(:, noiseSelect);
SRecAllSelect = SRecAll(:, driftSelect, noiseSelect);
% obtain the maximum of all reconstructions, so that the colormap is good
maxVecFun = @(x) max(x(:));  % max of a vector
arrayOfMaxCellFun = @(x) cellfun(maxVecFun,  x); % exact max of each cell element, form an maximum array (assume each cell element is an array)
maxI = max(arrayOfMaxCellFun(WUncalibAllSelect(:))); % 1.50, but maxI = max(arrayOfMaxCellFun(WUncalibAllSelect(:))) is 360
% but the maxI = 1.50 make the figure too dark, so still manually set it to 1
maxI = 1.0;

cMap = parula(2^16);  % convert grayscale [0, 1] to floating rgb with cMap
clim = [0, maxI]; 
%climSino = [0, max(arrayOfMaxCellFun([SDriftAllSelect(:); SNoDriftAllSelect(:)]))];
climSinoAll = {[0, max(SNoDriftAll{1,1}(:))];  [0, max(SNoDriftAll{2,1}(:))]}; % two sample have different maximum sinogram intensity
 
myFontSize = 32;

% generate L matrix without drift for sinogram
WSz = size(WRecAll{1}); [NTheta, NTau] = size(SDriftAll{1});
L = XTM_Tensor_XH(WSz, NTheta, NTau);
% only L matrix on or after 20190125 is normalized in simulation
LNormalizer = full(max(sum(L,2)));  s0 = SNoDriftAll{1,1}(:); s = L*WGTAll{1}(:); % Compute the normalizer for 1st scan of no drift so that maximum row sum is 1 rather than a too small number
if norm(s0-s*LNormalizer) < norm(s0-s); L = L/LNormalizer; end

mean(SNoDriftAll{1,1}(:))

for nSample = 1:NSample
    sampleName = sampleNameAll{nSample};   
    climSino = climSinoAll{nSample};    
    WGT = WGTAll{nSample};
    
    imgName = sprintf('%s_%s', retFile(5:17), sampleName); % retFile(5:17) is timeTag    
    imwrite(gray2rgb(WGT, cMap, clim), fullfile(retFolder, [imgName '_WGT.png']));
    
    for nMaxDrift = 1:numel(maxDriftSelect)
        maxDrift = maxDriftSelect(nMaxDrift);                
                   
        for nNoise = 1:numel(gaussianSTDAllSelect)
            gaussianSTD = gaussianSTDAllSelect(nNoise);
            %%% sinogram no drift
            sinoObserved = SNoDriftAllSelect{nSample, nNoise};
            imgName = sprintf('%s_%s_noise%.3f', retFile(5:17), sampleName, gaussianSTD); % retFile(5:17) is timeTag
            imgName = regexprep(imgName, '\.(\d+)', '_$1'); % replace .12 to _12 for latex include figure files                                
            imwrite(gray2rgb(sinoObserved, cMap, climSino), fullfile(retFolder, [imgName '_sinoObserved.png']));
%             figure(1); clf; imagesc(sinoObserved, climSino); colormap(cMap);axis image; xlabel('\tau'); ylabel('\theta'); set(gca,'xtick',[]); set(gca,'ytick',[]);
%             printfig(fullfile(retFolder, [imgName '_sinoObserved_withLabel.png']), [800, 400], 'fontSize', myFontSize, 'overwrite', 1);
            %%% results with drift and noise
            % drift = driftRecAllSelect{nSample, nMaxDrift, nNoise};
            if isUsedBestWtUncalib
                WUncalib = WUncalibBestWtAllSelect{nSample, nMaxDrift, nNoise}; % best weight
            else
                WUncalib = WUncalibAllSelect{nSample, nMaxDrift, nNoise}; % fixed weight
            end
            

            WRec = WRecAllSelect{nSample, nMaxDrift, nNoise};
            sinoObserved = SDriftAllSelect{nSample, nMaxDrift, nNoise} ;
            SRec = SRecAllSelect{nSample, nMaxDrift, nNoise}; 
        
            imgName = sprintf('%s_%s_maxDrift%.1f_noise%.3f', retFile(5:17), sampleName, maxDrift, gaussianSTD); % retFile(5:17) is timeTag
            imgName = regexprep(imgName, '\.(\d+)', '_$1'); % replace .12 to _12 for latex include figure files
            fprintf('%s:  psnr_rec=%.2fdB, psnr_baseline=%.2fdB.\n', imgName, psnr(WRec, WGT), psnr(WUncalib(:,:,end), WGT));
           
            imwrite(gray2rgb(WRec, cMap, clim), fullfile(retFolder, [imgName '_ours.png']));
            imwrite(gray2rgb(WUncalib(:,:,end), cMap, clim), fullfile(retFolder, [imgName '_baseline.png']));
%             figure(1);clf; imagesc(WRec, clim); axis image; %xlabel('x'); ylabel('y');            
%             printfig(fullfile(retFolder, [imgName '_oursOld.png']), [1000, 1000])
%             figure(1);clf; imagesc(WUncalib(:,:,end), clim); axis image; %xlabel('x'); ylabel('y');            
%             printfig(fullfile(retFolder, [imgName '_baselineOld.png']), [1000, 1000])
            
            % sinogram
            imwrite(gray2rgb(reshape(L*WRec(:), [NTheta, NTau]), cMap, climSino), fullfile(retFolder, [imgName '_Rec2SinoWithNoDriftModel_ours.png']));
            imwrite(gray2rgb(reshape(L*reshape(WUncalib(:,:,end), [], 1), [NTheta, NTau]), cMap, climSino), fullfile(retFolder, [imgName '_Rec2SinoWithNoDriftModel_baseline.png']));            
            
            imwrite(gray2rgb(sinoObserved, cMap, climSino), fullfile(retFolder, [imgName '_sinoObserved.png']));
            %figure(1); clf; imagesc(sinoObserved, climSino); colormap(cMap);axis image; xlabel('\tau'); ylabel('\theta');      
            %printfig(fullfile(retFolder, [imgName '_sinoObserved.png']), [800, 400], 'fontSize', myFontSize, 'overwrite', 1); %                                    
            
        end
    end
end
%% save colormap
% figure, imagesc(WGT); colormap(cMap); axis image; colorbar;
% printfig(fullfile(retFolder, 'colormap.png'), [1200, 1000], 'fontSize', 16, 'zeroWhiteMargin', 0, 'overwrite', 1); % need to crop by hand afterwards
%% Determine the maximum intensity for each sample's reconstruction

