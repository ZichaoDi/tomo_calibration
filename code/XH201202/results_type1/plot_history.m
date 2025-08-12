function plot_history(matFile)
% plot_history  –  Simple comparison plot for PSNR, misfit, lambda and images
% Old = blue dashed w/ circles, new = red solid w/ diamond.
% Convergence markers on PSNR plot are placed at iteration where
% the change in misfit term is ≤ 1e-3.
% Also plots W, final WRecs, and final WRecs_old.

% --- Load necessary variables from the MAT file
load(matFile, 'info', 'info_old', 'W', 'WRecs', 'WRecs_old');

% --- Plotting for PSNR, Misfit, and Lambda history
fields = {'psnr_history', 'misfit_history', 'lambda_history'};
yLabs  = {'PSNR(dB)', 'Mis-fit term', '\lambda'};
x      = 1:10;

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
        else
            convIter_old = numel(x); % fallback to last iteration
        end
        h3 = plot(x(convIter_old), yOld(convIter_old), 'bo', 'MarkerFaceColor', 'b', ...
                  'MarkerSize', 16, 'LineWidth', 4);
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        % --- New approach convergence point: Red diamond
        delta_new = abs(diff(info.misfit_history(1:numel(x))));
        idx_new = find(delta_new <= 1e-3, 1, 'first');
        if ~isempty(idx_new) && idx_new + 1 <= length(yNew)
            convIter_new = idx_new + 1;
        else
            convIter_new = numel(x); % fallback to last iteration
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

% --- Plotting for the three requested images without titles/axes ---

% 1st: Ground truth image
figure;
imagesc(W);
axis image off;
colorbar off;

% 2nd: Final reconstruction from the 'new' approach
figure;
imagesc(reshape(WRecs(:,end), 256, 256));
axis image off;
colorbar off;

% 3rd: Final reconstruction from the 'old' approach
figure;
imagesc(reshape(WRecs_old(:,end), 256, 256));
axis image off;
colorbar off;

end