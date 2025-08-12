function [tv, tvNormalized] = TVnorm3D(x, voxSz, TVChoice)
%TVNORM3D Summary of this function goes here
%   Detailed explanation goes here
%   hack on 2016.03.28, changed on 2016.04.05
%  2016.09.10_updated to include voxSz, because the gradient should be finite difference divide by step/voxel size
%           non-uniform voxSz has non-uniform step size of computing , however, the result with/without voxSz=[98,98,50] are
%           very similar(same) for Bacteria1_DoubleLobes_MFM_10Measures_frame=2
%  2016.11.07_include the isotropic version, and make the d1, d2, d3 size the same
%  2018.08.13: bring back TVnorm for 2D image

%% Naive 2D TV norm
if ndims(x) == 2 % 2D
    tv = sum(sum(sqrt(diffh(x).^2+diffv(x).^2)));
    return
end

%% 3D TV norm

if ~exist('TVChoice', 'var') || isempty(TVChoice)
    TVChoice = 'isotropic';
end

sz = size(x);
% 
% forward difference, problem: object of size 5x5x5 has TVnorm of 100, but L1norm of 125, close
d1 = x([2:sz(1), sz(1)], :, :) - x;
d2 = x(:, [2:sz(2), sz(2)], :) - x;
d3 = x(:, :, [2:sz(3), sz(3)]) - x;

if ~exist('voxSz', 'var') || isempty(voxSz) || all(voxSz == voxSz(1)) %all(voxSz==1) 
    % we have no information on step/voxel size, assume it is uniform.
    if exist('TVChoice', 'var') && strcmpi(TVChoice, 'anisotropic')
        tv = sum(abs(d1(:)) + abs(d2(:)) + abs(d3(:)));
    else
        % default istropic
        tv = sum(sqrt(d1(:).^2 + d2(:).^2 + d3(:).^2));
    end
else
    stepSz = voxSz/max(voxSz); % normlized stepSize
    if exist('TVChoice', 'var') && strcmpi(TVChoice, 'anisotropic')    
        tv = sum(abs(d1(:))/stepSz(1) + abs(d2(:))/stepSz(2) + abs(d3(:))/stepSz(3));
    else
         % default istropic
        tv = sum(sqrt((d1(:)/stepSz(1)).^2 + (d2(:)/stepSz(2)).^2 + (d3(:)/stepSz(3)).^2));
    end
end

% normalize the sparsity with maximum absolute value and the number of elements
if nargout == 2
    %normalizer = max(abs(d1(:))/stepSz(1) +  abs(d2(:))/stepSz(2) + abs(d3(:))/stepSz(3)) * numel(d1);
    normalizer = (max(abs(d1(:))/stepSz(1)) +  max(abs(d2(:))/stepSz(2)) + max(abs(d3(:))/stepSz(3))) * numel(d1);
    tvNormalized = tv / normalizer;
end

%ret = sum(abs(d1(:))) + sum(abs(d2(:))) ;

% ret=0;
% for iz=1:size(x,3)
%     %ret = ret + TVnorm(x(:,:,iz));
%     x2D = x(:,:,iz);
%     ret = ret + sum(sum(sqrt(diffh(x2D).^2+diffv(x2D).^2)));
% end
end

% OLD version of handle boundary
%d1 = x(2:end, :, :) - x(1:end-1, :, :);
%d2 = x(:, 2:end, :) - x(:, 1:end-1, :);
%d3 = x(:, :, 2:end) - x(:, :, 1:end-1);

