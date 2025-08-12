function metrics = plot_history(matFile)
% plot_history  –  Plots history, returns key performance metrics.
%
% Plots comparison for PSNR, misfit, lambda and converged images.
% Returns a struct with PSNR/SSIM values at the point of convergence.

% --- Load necessary variables from the MAT file
load(matFile, 'info', 'info_old', 'W', 'WRecs', 'WRecs_old', 'PSNR_baseline', 'SSIM_baseline');

% --- Plotting for PSNR, Misfit, and Lambda history
fields = {'psnr_history', 'misfit_history', 'lambda_history'};
yLabs  = {'PSNR(dB)', 'Mis-fit term', '\lambda'};
x      = 1:10;

% Initialize convergence iterators before the loop to ensure they are in scope later
convIter_old = numel(x);
convIter_new = numel(x);

for k = 1:3
    yOld = info_old.(fields{k})(1:numel(x));
    yNew = info.(fields{k})(1:numel(x));
    
    figure; hold on;
    h1 = plot(x, yOld, 'bo--', 'LineWidth', 4, ...
              'MarkerSize', 7, 'MarkerFaceColor', 'w');
    h2 = plot(x, yNew, 'r-',   'LineWidth', 4);
    
    % --- Place convergence markers on PSNR plot ONLY (k==1)
    if k == 1
        % --- Old approach convergence point: Solid blue circle
        delta_old = abs(diff(info_old.misfit_history(1:numel(x))));
        idx_old = find(delta_old <= 1e-3, 1, 'first');
        if ~isempty(idx_old) && idx_old + 1 <= length(yOld)
            convIter_old = idx_old + 1;
        end
        h3 = plot(x(convIter_old), yOld(convIter_old), 'bo', 'MarkerFaceColor', 'b', ...
                  'MarkerSize', 16, 'LineWidth', 4);
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        % --- New approach convergence point: Red diamond
        delta_new = abs(diff(info.misfit_history(1:numel(x))));
        idx_new = find(delta_new <= 1e-3, 1, 'first');
        if ~isempty(idx_new) && idx_new + 1 <= length(yNew)
            convIter_new = idx_new + 1;
        end
        h4 = plot(x(convIter_new), yNew(convIter_new), 'rd', ...
                  'MarkerFaceColor', 'r', 'MarkerSize', 16, 'LineWidth', 4);
        set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    
    xlabel('Iteration k'); ylabel(yLabs{k});
    legend([h1 h2], {'old approach','new approach'}, ...
           'Location','southoutside','Orientation','horizontal');
    xlim([1 10]); set(gca,'FontSize',16); grid off;
end

% --- Plotting for the three requested images at convergence ---
% 1st: Ground truth image
figure;
imagesc(W);
axis image off;
colorbar off;

% 2nd: Reconstruction from 'new' approach at convergence
figure;
imagesc(reshape(WRecs(:,convIter_new), 256, 256));
axis image off;
colorbar off;

% 3rd: Reconstruction from 'old' approach at convergence
figure;
imagesc(reshape(WRecs_old(:,convIter_old), 256, 256));
axis image off;
colorbar off;

% --- Prepare the output struct with the requested metrics ---
metrics.psnr_new = info.psnr_history(convIter_new);
metrics.psnr_old = info_old.psnr_history(convIter_old);
metrics.psnr_baseline = PSNR_baseline;

% Check if ssim_history exists to avoid errors
if isfield(info, 'ssim_history')
    metrics.ssim_new = info.ssim_history(convIter_new);
else
    metrics.ssim_new = NaN; % Use NaN as a placeholder
end

if isfield(info_old, 'ssim_history')
    metrics.ssim_old = info_old.ssim_history(convIter_old);
else
    metrics.ssim_old = NaN; % Use NaN as a placeholder
end

metrics.ssim_baseline = SSIM_baseline;

end