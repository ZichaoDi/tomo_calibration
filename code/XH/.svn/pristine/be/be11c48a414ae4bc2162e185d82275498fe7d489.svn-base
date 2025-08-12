function W = createObject(sampleName, WSz)
%Create a sampleName object W 
%Inputs:
% sampleName:   the name that specifies a sample, default 'Phantom", not case-sensitive
% WSz:          size of W, 2D array [Ny, Nx].
%Output: 
% W:            object volume of size [Ny, Nx]
% 
% 2018.08.17:   XH refactored from script to function
% 2018.08.20:   XH refactored to eliminate duplcated inputs WSz and WSz(1), moved x, y and omega out of this function
% 2018.11.28:   XH cleaned, and rename to createObject for a general



%%%=========== assign weight matrix for each element in each pixel
% Z: Atomic Numbers of Existing Elements in the object
switch lower(sampleName)
    case {'golosio', 'checkboard'}
        Z = 19;
    case {'phantom', 'brain'}
        Z = 19;
    case 'circle'
        Z = 14;
    otherwise
        error('Unknown Sample Specification: %s!', sampleName)
end


%---------------------------
NumElement = length(Z);
W = zeros(WSz(1), WSz(2), NumElement);
switch lower(sampleName)
    case {'phantom'}
        W = CreateElement_XH(WSz, NumElement);
    case 'brain'        
        W0 = imread('MRI_T2_Brain_axial_image.png'); % original source converted to .png:  https://commons.wikimedia.org/wiki/File:MRI_T2_Brain_axial_image.jpg
        W = resizeWithSameAspectRatio(W0, WSz);
    case 'golosio'
        CreateCircle; % XH: may need rewrite if it is used
    case 'circle'
        [X,Y]=meshgrid(1:WSz(2), 1:WSz(1)); % (1:XMax, 1:YMax)
        center=[WSz(2)/2, WSz(1)/2]; % (center_x, center_y)
        r=20;
        pix = (X-center(1)).^2+(Y-center(2)).^2 <= r^2;%& (X-center(1)).^2+(Y-center(2)).^2>=(r-5).^2; %% circle
        W(pix)=10;
    case 'checkboard'
        W = kron(invhilb(WSz(1)/10)<0, ones(10,10));
        W = repmat(W,[1 1 NumElement]);        
    case lower('fakeRod')
        CreateRod;  % XH: may need rewrite if it is used
    otherwise
        error('Unknown Sample Specification: %s!', sampleName)
end

assert(all(size(W)==WSz), 'The output size of W do not match the user specified size by WSz!');

% Normalize pixel value to 0 and 1
W(W<0) = 0;
W = W/max(W(:));
%W = W*40;  %TODO: DEBUG % 2018.09.12_Wendy added *40 for a bug of the solver tnbc(), 2018.12.21_make the maximum intensity 1,
% temp fix the problem by scale the data inside the solver by *40.

end



function W = resizeWithSameAspectRatio(W0, WSz)
% resize that preserve the image aspect ratio
% 1. If the size of original image W0 is different from user specified size WSz, we resize W0 to WSz using bicubic interpolation
% 2. However, if the aspect ratio of WOrig is different from WSz, we pad W0 with repeated boundry to achieve the same ratio before resize

W0Sz = size(W0);

% (1) Check aspect ration, if not much, pad W0
sz1 = round(WSz(1)/WSz(2) * W0Sz(2));
sz2 = round(WSz(2)/WSz(1) * W0Sz(1));
if (W0Sz(1) ~= sz1) ||  (W0Sz(2) ~= sz2) % Compare in both direction to make sure, to avoid strange round effect
    % Aspect ratio different, we pad W0 to new size
    % We can't have both sz1>W0Sz(1) and sz2>W0Sz(2), and can't have both sz1<W0Sz(1) and sz2<W0Sz(2). 
    % Hence W0NewSz will be either [W0Sz(1), sz2] or [sz1, W0Sz(2)]
    W0NewSz = max(W0Sz, [sz1, sz2]);  % same as W0NewSz = max([W0Sz(1), sz2], [sz1, W0Sz(2)]); % elementwise maximum
    padSz = W0NewSz - W0Sz;  % only one of padSz(1) and padSz(2) is non-zero
    padSzPre = round(padSz/2);
    padSzPost = padSz - padSzPre;
    W0Pad = padarray(W0, padSzPre, 'replicate','pre');
    W0Pad = padarray(W0Pad, padSzPost, 'replicate','post');
end

% (2) resize W0 to size of WSz
W = imresize(double(W0Pad), WSz, 'bicubic'); % convert to double for any 


end