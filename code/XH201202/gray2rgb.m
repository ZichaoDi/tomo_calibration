function Irgb = gray2rgb(Igray, cMap, clim)
%Irgb = gray2rgb(Igray, cMap, clim)
%Convert a grayscale image to a Irgb color image, using a given colormap cMap (default is parula if not specified). 
% If clim = [cMin, cMax] is given, then intensity <= cMin will be mapped to first element of colormap, and intensity >= cMax
% will be mapped to the last element of colormap.  By default clim is [0-255] for uint8, [0-65535] for uint16 and [0, 1] for float.
%
% The purpose of this function is to make details of a grayscale image more visible, but convert it to a pseudo color rgb image
%
% Usage:
%   Irgb = gray2rgb(Igray)
%   Irgb = gray2rgb(Igray, 'parula');
%   Irgb = gray2rgb(Igray, 'jet');
%   Irgb = gray2rgb(Igray, 'jet', [0, 0.8]);
%
% Inputs: 
%   Igray: the grayscale image (2D for now), which could be uint8, uint16 or float format. Igray could also be a file name where
%   the image content is there
%   cMap: color map, 'parula'(default), 'jet', n*3 array
%   clim: [cMin, cMax] two element array, similar to clim in imagesc.
%       if clim is not specified, we assume that we map 
%
% References:
%   ind2rgb
%
%  Xiang Huang: xianghuang@gmail.com
%   2019.01.23 created 
% 
if nargin==0; testMe(); end

%% Parse Inputs
% 1. Parse Igray
if ischar(Igray)
    Igray = imread(Igray);
end
% Note: used ndims not ismatrix, though ismatrix is recommended by matlab, because it is confusing as ismatrix(rand(2,3,4)) = false for 3D matrix.
assert(ndims(Igray)==2 && any(ismember({'uint8', 'uint16', 'single', 'double'}, class(double(Igray)))),...
    'input must be a 2D image of uint8/uint16 interger, or single/double floating point format');  

% 2. Parse cMap
if ~exist('cMap', 'var') || isempty(cMap)
    cMap = 'parula';  % default parula
end
if ischar(cMap)
    switch class(Igray)
        case 'uint8'
            NColors = 2^8;
        case {'uint16', 'single', 'double'}
            NColors = 2^16;        
        otherwise
            error('unknown numerical class of gray image');            
    end
    cMap = feval(cMap, NColors);  % default 256 color (matlab default is 64 which is bad)
else
    assert(size(cMap,1)==2^8 || size(cMap,1)==2^16, 'currently user specified colormap must has 2^8 or 2^16 rgb triplets.');
end
[NColors, NChannel] = size(cMap);
assert(isnumeric(cMap) && NChannel==3, 'cMap must be N*3 color array or a valid colormap string such as parula, jet'); 

% 3. Parse clim
if ~exist('clim', 'var') || isempty(clim)
    switch class(Igray)
        case 'uint8'
            clim = [0, 2^8-1];
        case 'uint16'
            clim = [0, 2^16-1];
        case {'single', 'double'}
            clim = [0, 1];            
        otherwise
            error('unknown numerical class of gray image');            
    end
end


%% Map Igray intensity from clim = [cMin, cMax] to index of [1, NColors]
% 1. convert to double and clamp to [cMin, cMax]
Igray = double(Igray);
cMin = clim(1);
cMax = clim(2);
Igray(Igray<cMin) = cMin;
Igray(Igray>cMax) = cMax;
% 2. Map from [cMin, cMax] intensity image to [1, NColors] integer index image
% For example, map [0, 255] uint8 grayscale iamge to [1, 256] index image, i.e., Iidx = Igray + 1;
Iidx = round((NColors-1)/(cMax-cMin) * (Igray-cMin) + 1);

% 3. Index image to rgb image
Irgb = ind2rgb(Iidx, cMap);


end


function testMe()
%% uint8
Igray = imread('cameraman.pgm');
cMap = 'parula';
Irgb = gray2rgb(Igray, cMap);
figure, clf; imagesc(Igray); colormap(cMap); axis image; colorbar;
figure, clf; imshow(Irgb); axis image; colorbar;

%% floating
Igray = imread('cameraman.pgm'); Igray = double(Igray)/255;
cMap = 'parula';
Irgb = gray2rgb(Igray, cMap);
figure, clf; imagesc(Igray); colormap(cMap); axis image; colorbar;
figure, clf; imshow(Irgb); axis image; colorbar;


%% climed
Igray = imread('cameraman.pgm'); Igray = double(Igray)/255;
cMap = 'parula';
Irgb = gray2rgb(Igray, cMap, [0, 0.5]);
figure, clf; imagesc(Igray); colormap(cMap); axis image; colorbar;
figure, clf; imshow(Irgb); axis image; % colorbar; %colorbar of rgb makes no sense?
end