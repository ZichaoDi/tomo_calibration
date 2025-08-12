%% 1st figure
% Row 1: object:   Groundtruth, Recontruction without Calibration, Reconstruction with Calibration (ours), (Reconstruction with
% Know Drift, Previous result if any)
% Row 2: sinogram: Corresponding sinagrom of first row

%WAllSave = {WGT, WWendy(:,:,3), X}; 
%WAllSaveNames = {'obj_GT', 'obj_rec_uncalibrated', 'obj_rec_calibrated'};
% todo change WWendy to WXH

timeTag = '20181113';
WAllSave = {WGT, X,  WWendy(:,:,3), WWendy(:,:,1), WWendy(:,:,2), WXH(:,:,3), WXH(:,:,1), WXH(:,:,2)};
WAllSaveNames = {'obj_GT', 'obj_rec_calibrated_10Iter', ...
                 'obj_rec_uncalibrated_wendy', 'obj_rec_no_drift_wendy', 'obj_rec_knownCalibration_wendy', ...
                 'obj_rec_uncalibrated_twist', 'obj_rec_no_drift_twist', 'obj_rec_knownCalibration_twist'};

for n = 1:numel(WAllSaveNames)
    psnr_= difference(WAllSave{n}, WAllSave{1}); 
    WAllSaveNames{n} = sprintf('%s_%s_%.2fdB', WAllSaveNames{n}, timeTag, psnr_);
end
WAllSaveNames = cellfun(@(s) regexprep(s, '\.(\d+)dB', '_$1dB'), WAllSaveNames,'UniformOutput',false); % replace .12dB to _12dB for latex include figure files

SAllSave = {SAll(:,:,2), SRec, ...
            reshape(LAll{1}*reshape(WWendy(:,:,3), [], 1), NTheta, NTau), ...
            reshape(LAll{1}*reshape(WWendy(:,:,1), [], 1), NTheta, NTau), ...
            reshape(LAll{2}*reshape(WWendy(:,:,2), [], 1), NTheta, NTau), ...
            reshape(LAll{1}*reshape(WXH(:,:,3), [], 1), NTheta, NTau), ...
            reshape(LAll{1}*reshape(WXH(:,:,1), [], 1), NTheta, NTau), ...
            reshape(LAll{2}*reshape(WXH(:,:,2), [], 1), NTheta, NTau), ...
            SAll(:,:,1)}; 
SAllSaveNames = {'sino_GT', 'sino_calibrated_10Iter', ...
                 'sino_uncalibrated_wendy', 'sino_no_drift_wendy', 'sino_knownCalibration_wendy', ...
                 'sino_uncalibrated_twist', 'sino_no_drift_twist', 'sino_knownCalibration_twist', ...
                 'sino_no_drift_GT'};

for n = 1:numel(SAllSaveNames)
    SAllSaveNames{n} = sprintf('%s_%s_%1.2fdB', SAllSaveNames{n}, timeTag, difference(SAllSave{n}, SAllSave{1}));
end
SAllSaveNames = cellfun(@(s) regexprep(s, '\.(\d+)dB', '_$1dB'), SAllSaveNames,'UniformOutput',false); % replace .12dB to _12dB for latex include figure files

retFolder = fullfile('..', '..', 'writing', 'driftedTomography', 'figs'); %'xx\Dropbox\Tomography\writing\driftedTomography\figs';


save(fullfile(retFolder, ['ret_', timeTag]), 'WAllSave', 'SAllSave', 'WAllSaveNames', 'SAllSaveNames', 'retFolder', 'timeTag');

maxVecFun = @(x) max(x(:));  % max of a vector
arrayOfMaxCellFun = @(x) cellfun(maxVecFun,  x); % exact max of each cell element, form an maximum array (assume each cell element is an array)


%% Save figures for row 1: object
maxI = max(arrayOfMaxCellFun(WAllSave(1:3))); % 92.4078, same as maxI = max(arrayOfMaxCellFun(WAllSave(1:5)))

figure(999);clf; 
for n = 1:numel(WAllSave)
    imagesc(WAllSave{n}, [0, maxI]); axis image; %xlabel('x'); ylabel('y');
    %xlabel('x'); ylabel('y');
    printfig(fullfile(retFolder, [WAllSaveNames{n} '.png']), [1000, 1000])
end
colorbar; printfig(fullfile(retFolder, ['obj_colorbar_', timeTag, '.png']), [1400, 1000])


%% Save figures for row 2: sinogram
maxI = max(arrayOfMaxCellFun(SAllSave(1:3))); % 0.2030, similar to 0.2040 maxI = max(arrayOfMaxCellFun(SAllSave(1:5)))

figure(999);clf; 
for n = 1:numel(SAllSave)
    imagesc(SAllSave{n}, [0, maxI]); axis image; xlabel('\tau'); ylabel('\theta');
    %xlabel('x'); ylabel('y');
    printfig(fullfile(retFolder, [SAllSaveNames{n} '.png']), [1000, 1000])
end
colorbar; printfig(fullfile(retFolder, ['sino_colorbar_', timeTag, '.png']), [1400, 1000])
%%

figure(991); clf, multAxes(@imagesc, WAllSave); multAxes(@axis, 'image'); linkAxesXYZLimColorView(); 
printfig(fullfile(retFolder, ['obj_zzz_all_', timeTag, '.png']), [2560, 1600]);
figure(992); clf; multAxes(@imagesc, SAllSave); multAxes(@axis, 'image'); linkAxesXYZLimColorView();
printfig(fullfile(retFolder, ['sino_zzz_all_', timeTag, '.png']), [2560, 1600]);