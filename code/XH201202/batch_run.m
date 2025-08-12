%% batch_run_tomo.m  –  fire-and-forget batch runner
% --------------------------------------------------
%  EDIT THIS LINE if your main file has a different name
% --------------------------------------------------
MAIN_FILE = 'mainScript1.m';     % <-- change to your script or function
% -------- parameter grids -----------------------------------------------
sampleNames   = {'Brain','Phantom'};
gaussSTDs     = [0,0.01,0.02];
maxDrifts     = [1 2 3 5];
lambda0s      = 1e-6;
% -------- output folder --------------------------------------------------
outDir = fullfile(pwd,'results_type2');
if ~exist(outDir,'dir'); mkdir(outDir); end
% -------- loops ----------------------------------------------------------
for iS = 1:numel(sampleNames)
for iN = 1:numel(gaussSTDs)
for iD = 1:numel(maxDrifts)
for iL = 1:numel(lambda0s)
    % House-keeping
    clearvars -except MAIN_FILE outDir                               ...
                         sampleNames gaussSTDs maxDrifts lambda0s    ...
                         iS iN iD iL
    close all force
    % Variables the main script/function will read
    batchMode   = true;                         
    sampleName  = sampleNames{iS};             
    gaussianSTD = gaussSTDs(iN);                
    maxDrift    = maxDrifts(iD);                
    lambda0     = lambda0s(iL);                 
    fprintf('\n>>> %-7s | σ=%.3f | drift=%d | λ₀=%g\n', ...
            sampleName, gaussianSTD, maxDrift, lambda0);
    % ---------------------------------------------------------------------
    %  CALL YOUR RECONSTRUCTION CODE
    % ---------------------------------------------------------------------
    % A) if MAIN_FILE is a ***script***:
    run(MAIN_FILE);
    % B) if MAIN_FILE is a ***function***, comment out A) and use:
    % [WRecs,SRecs,drift,driftAll,info, ...
    %  WRecs_old,SRecs_old,drift_old,driftAll_old,info_old, ...
    %  W, PSNR_baseline, SSIM_baseline] = ... % Added W, PSNR_baseline, SSIM_baseline here
    %     mainScript1(batchMode,sampleName,gaussianSTD,maxDrift,lambda0);
    % ---------------------------------------------------------------------
    %  ALWAYS save the exact ten matrices/structs requested
    % ---------------------------------------------------------------------
    % Updated: Added 'W', 'PSNR_baseline', 'SSIM_baseline' to varsWanted
    varsWanted = { ...
        'WRecs','SRecs','drift','driftAll','info', ...
        'WRecs_old','SRecs_old','drift_old','driftAll_old','info_old', ...
        'W', 'PSNR_baseline', 'SSIM_baseline'}; % New variables added
    % Ensure every name exists; create empty placeholder if not
    for k = 1:numel(varsWanted)
        if exist(varsWanted{k},'var') ~= 1
            eval([varsWanted{k} ' = [];']);   
        end
    end
    % File name
    outFile = fullfile(outDir, ...
        sprintf('result_%s_sigma%.3g_drift%d_lambda%.0e.mat', ...
                sampleName, gaussianSTD, maxDrift, lambda0));
    % Save
    % Updated: Included new variables in the save command
    save(outFile, varsWanted{:}, ...
         'sampleName','gaussianSTD','maxDrift','lambda0','-v7.3');
    fprintf('    saved → %s\n', outFile);
end
end
end
end

