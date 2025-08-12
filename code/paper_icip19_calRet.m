% script of running result with multiple settings
% NOTE: noise really should be Poisson plus Guassian on sinogram
% NOTE: we currently run once for one maxDrift and one noise, strictly may need to run multiple times for each maxDrift/noise
% and compute average
% NOTE: we used 1e-8 regularizer weight for all noise level, 1e-8 is good for noiseless case, the optimal regWt should be larger for noisy case
% NOTE: assume maxDrift is known, but in reality may not know?
% 20190124:19pm swapped the order of nNoise, nMaxDrift to nMaxDrift, nNoise

opengl('save','software'); %hack fix for crashed matlab due to the graphics driver; https://www.mathworks.com/matlabcentral/answers/157894-resolving-low-level-graphics-issues-in-matlab
batchMode = true; % turn on batch mode, so that sampleName, maxDrift and noiseStd is set here rather than in mainXH_XH
sampleNameAll = {'Brain', 'Phantom'}; % choose from {'Phantom', 'Brain', ''Golosio', 'circle', 'checkboard', 'fakeRod'}; not case-sensitive
gaussianSTDAll = [0, 0.01, 0.02]; % [0, 0.01, 0.02, 0.05];  % noise level, 1st run wrong! noise should be add to sinogram not W

maxDriftAll = [1, 2, 3, 5];  % % maximum error; need to add 0.5 for floating
maxMaxDriftAll = max(maxDriftAll);
NSample = numel(sampleNameAll);
NNoise = numel(gaussianSTDAll);
NMaxDrift = numel(maxDriftAll);

%%
WGTAll = cell(NSample, 1);   % may change to WGTAll = cell(NSample, 1); Ground truth object, with different noise level added, should be indepdnet on noise leve and drift
SNoDriftAll = cell(NSample, NNoise); % Sinogram without drift, depend on noise, independent of drift
cell0 = cell(NSample, NMaxDrift, NNoise); 
SDriftAll = cell0; % Sinogram with drift, depends on both noise and dirft
driftGTAll = cell0; % ground truth drift, for each noise and maxDrift

% Reconstruction 
driftRecAll = cell0;
WUncalibAll  = cell0; % uncalibrated reconstruction: known-no-drift, known-with-drift, unkown-with-drift
WRecAll = cell0; % calibrated reconstruction
SRecAll = cell0;
regWtCalibAll = cell0;

% 2019.01.24_uncalibrated reconstruction with tries of different weights
WUncalibBestWtAll = cell0;
regWtBestUncalibAll = cell0;
psnrAllWtUncalibAll = cell0;
regWtAllUncalibAll = cell0;

for nSample = 1:NSample
    sampleName = sampleNameAll{nSample};
    for nMaxDrift = 1:NMaxDrift
        maxDrift = maxDriftAll(nMaxDrift);        
        for nNoise = 1:NNoise
            gaussianSTD = gaussianSTDAll(nNoise);
            % compute the result for each setting, by call the mainScript
            mainScript
            WGTAll{nSample} = WGT;
            SNoDriftAll{nSample, nNoise} = SAll(:,:,1);
            SDriftAll{nSample, nMaxDrift, nNoise} = SAll(:,:,2);
            driftGTAll{nSample, nMaxDrift, nNoise} = driftGT;
            
            % Reconstruction
            driftRecAll{nSample, nMaxDrift, nNoise} = drift;
            WUncalibAll{nSample, nMaxDrift, nNoise} = WUncalib;
            WRecAll{nSample, nMaxDrift, nNoise} = X;
            SRecAll{nSample, nMaxDrift, nNoise} = SRec;
            regWtCalibAll{nSample, nMaxDrift, nNoise} = regWt; 
            
            % 2019.01.24_uncalibrated reconstruction with tries of different weights
            WUncalibBestWtAll{nSample, nMaxDrift, nNoise} = WUncalibBestWt;
            regWtBestUncalibAll{nSample, nMaxDrift, nNoise} = regWtBestUncalib;
            psnrAllWtUncalibAll{nSample, nMaxDrift, nNoise} = psnrAllWtUncalib;
            regWtAllUncalibAll{nSample, nMaxDrift, nNoise} = regWtAllUncalib;
            
            close all;
        end
    end
end
curTime = datestr(now, 'yyyymmdd_HHMM');
save(['ret_' curTime '.mat'], ...
    'sampleNameAll', 'gaussianSTDAll','maxDriftAll', ...
    'WGTAll', 'SNoDriftAll', 'SDriftAll', 'driftGTAll', 'driftRecAll',  'WUncalibAll', 'WRecAll', 'SRecAll', ...
    'WUncalibBestWtAll', 'regWtBestUncalibAll', 'psnrAllWtUncalibAll', 'regWtAllUncalibAll', ...
    'regWt', 'regWtCalibAll', ...
    'maxIter', 'scaleBegin', 'scaleEnd'), 
%50 minutes for NSample*NNoise*NMaxDrift = 2*3*4 = 24, 2minutes each reconstruction
