% plot_and_save_all.m
% Generate plots, collect key metrics, and save results.

% ------------------ Parameter Grid --------------------------------------
sampleNames = {'Phantom', 'Brain'};
gaussSTDs   = [0, 0.01, 0.02];
maxDrifts   = [1, 2, 3, 5];
lambda0s    = 1e-3;

% ------------------ Output Folder ---------------------------------------
outDir = fullfile(pwd, 'figures');

if ~exist(outDir, 'dir'); mkdir(outDir); end

% ------------------ Initialize Results Collector ------------------------
results_table = table();

% ------------------ Loop over combinations ------------------------------
fprintf('Starting processing...\n');
for iS = 1:numel(sampleNames)
    for iG = 1:numel(gaussSTDs)
        for iD = 1:numel(maxDrifts)
            for iL = 1:numel(lambda0s)
                sampleName  = sampleNames{iS};
                sigma       = gaussSTDs(iG);
                maxDrift    = maxDrifts(iD);
                lambda0     = lambda0s(iL);

                matFile = sprintf('result_%s_sigma%.3g_drift%d_lambda%.0e.mat', ...
                                  sampleName, sigma, maxDrift, lambda0);

                if ~isfile(matFile)
                    fprintf('Skipped: %s (not found)\n', matFile);
                    continue;
                end

                fprintf('Processing: %s\n', matFile);
                % try
                    % Call plotting function and get metrics back
                    metrics = plot_history(matFile);

                    % Add parameter info to the metrics struct
                    current_run = table({sampleName}, sigma, maxDrift, 'VariableNames', ["sampleName", "sigma", "maxDrift"]);
                    metrics_table = struct2table(metrics);
                    
                    % Append the full row to our main results table
                    results_table = [results_table; [current_run, metrics_table]];

                % catch err
                %     fprintf('⚠️ Error processing %s: %s\n', matFile, err.message);
                %     close all; continue;
                % end

                % --- Save each open figure ---
                figs = findobj('Type', 'figure');
                for f = figs'
                    figName = sprintf('%s_sigma%.3g_drift%d_lambda%.0e_fig%d.svg', ...
                        sampleName, sigma, maxDrift, lambda0, f.Number);
                    outPath = fullfile(outDir, figName);
                    set(f, 'PaperPositionMode', 'auto');
                    print(f, outPath, '-dsvg');
                end
                close all;
            end
        end
    end
end

% ------------------ Save the final collected metrics --------------------
if ~isempty(results_table)
    % Save as a CSV file for easy viewing (e.g., in Excel)
    csv_path = fullfile(outDir, 'summary_metrics.csv');
    writetable(results_table, csv_path);
    fprintf('\n✅ Metrics successfully saved to: %s\n', csv_path);

    % Save as a .mat file for easy loading back into MATLAB
    mat_results_path = fullfile(outDir, 'summary_metrics.mat');
    save(mat_results_path, 'results_table');
    fprintf('✅ Metrics successfully saved to: %s\n', mat_results_path);
else
    fprintf('\nNo files were processed. No summary file was created.\n');
end