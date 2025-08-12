function plot_history(matFile)
% plot_history  –  Simple comparison plot for PSNR, misfit, lambda
% Old = blue dashed w/ circles, new = red solid w/ diamond.
% Diamond on PSNR plot is placed at iteration where convergence occurs
% (Δ misfit ≤ 1e-3). Only two legend entries shown.

load(matFile,'info','info_old');

fields = {'psnr_history','misfit_history','lambda_history'};
yLabs  = {'PSNR(dB)','Mis-fit term','\lambda'};
x      = 1:10;

for k = 1:3
    yOld = info_old.(fields{k})(1:numel(x));
    yNew = info.(fields{k})(1:numel(x));

    figure; hold on

    h1 = plot(x, yOld, 'bo--', 'LineWidth', 4, ...
              'MarkerSize', 7, 'MarkerFaceColor', 'w');
    h2 = plot(x, yNew, 'r-',   'LineWidth', 4);

    % --- Old final point: solid blue circle
    h3 = plot(x(end), yOld(end), 'bo', 'MarkerFaceColor', 'b', ...
              'MarkerSize', 16, 'LineWidth', 4);

    % --- Red diamond placement: depends on field
    if k == 1  % PSNR plot
        delta = abs(diff(info.misfit_history(1:numel(x))));
        delta
        idx = find(delta <= 1e-3, 1, 'first');  % first convergence point
        if ~isempty(idx) && idx+1 <= length(yNew)
            convIter = idx + 1;
        else
            convIter = numel(x);  % fallback to last
        end
        h4 = plot(x(convIter), yNew(convIter), 'rd', ...
                  'MarkerFaceColor', 'r', 'MarkerSize', 16, 'LineWidth', 4);
    end

    % Remove extras from legend
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    xlabel('Iteration k'); ylabel(yLabs{k});
    legend([h1 h2], {'old approach','new approach'}, ...
           'Location','southoutside','Orientation','horizontal');
    xlim([1 10]); set(gca,'FontSize',16); grid off
end
end
