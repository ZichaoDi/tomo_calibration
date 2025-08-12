% plot_and_save_all.m
% Generate PSNR/Misfit/Lambda plots and image plots from .mat files
% and save them as high-res vector graphics (SVG).

% Figure 1: PSNR history (psnr_history)
% Figure 2: Mis-fit term history (misfit_history)
% Figure 3: Lambda history ( lambda _history )
% Figure 4: Ground truth image ( imagesc (W) )
% Figure 5: Final 'new' reconstruction ( imagesc(reshapeWRecs(:,end) ,256,256)) )
% Figure 6: Final 'old' reconstruction ( imagesc(reshape (WRecs_old(:, end), 256, 256)) )

% ------------------ Parameter Grid --------------------------------------
sampleNames = {'Phantom', 'Brain'};
gaussSTDs   = [0, 0.01, 0.02];
maxDrifts   = [1, 3, 5];
lambda0s    = 1e-3;

% ------------------ Output Folder ---------------------------------------
outDir = fullfile(pwd, 'figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end

% ------------------ Loop over combinations ------------------------------
for iS = 1:numel(sampleNames)
    for iG = 1:numel(gaussSTDs)
        for iD = 1:numel(maxDrifts)
            for iL = 1:numel(lambda0s)
                sampleName  = sampleNames{iS};
                sigma       = gaussSTDs(iG);
                maxDrift    = maxDrifts(iD);
                lambda0     = lambda0s(iL);

                % Build MAT file name
                matFile = sprintf('result_%s_sigma%.3g_drift%d_lambda%.0e.mat', ...
                                  sampleName, sigma, maxDrift, lambda0);

                if ~isfile(matFile)
                    fprintf('Skipped: %s (not found)\n', matFile);
                    continue;
                end

                fprintf('Processing: %s\n', matFile);
                try
                    plot_history(matFile);  % Call the updated plotting function
                catch err
                    fprintf('⚠️ Error plotting %s: %s\n', matFile, err.message);
                    close all; continue;
                end

                % Save each open figure (now includes 3 line plots + 3 images)
                figs = findobj('Type', 'figure');
                for f = figs'
                    figName = sprintf('%s_sigma%.3g_drift%d_lambda%.0e_fig%d.svg', ...
                        sampleName, sigma, maxDrift, lambda0, f.Number);
                    outPath = fullfile(outDir, figName);
                    set(f, 'PaperPositionMode', 'auto');
                    print(f, outPath, '-dsvg'); % Save as SVG
                end
                close all; % Close all figures before the next iteration
            end
        end
    end
end