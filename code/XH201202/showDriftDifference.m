function showDriftDifference(drift, driftGT, mask)


if ~exist('mask', 'var') || isempty(mask)
    mask = true(size(drift)); 
end


driftGT_masked = driftGT; driftGT_masked(~mask) = 0; 
drift_masked = drift; drift_masked(~mask) = 0;
multAxes(@imagesc, {driftGT_masked, drift_masked, drift_masked-driftGT_masked}); multAxes(@title, {'drift-GT', 'drift-Rec', sprintf('psnr=%.2fdB', difference(drift_masked, driftGT_masked))}); linkAxesXYZLimColorView(); multAxes(@colorbar);
subplot(224); histogram(drift_masked-driftGT_masked, 50); title('hist of error')

% % old 1D comparison
% dGT=driftGT(mask); d=drift(mask); err_mean =  mean(abs(d(:)-dGT(:)));
% figure, plot(dGT, '-rx'); hold on; plot(d, '--bo'); title(sprintf('mean error =%.2e',err_mean)); legend('GT', 'calc'); grid on; ylim([-1.1, 1.1]);
% axes('Position',[.7 .7 .2 .2]); box on; histogram(d-dGT, 50); title('hist of error')


end