function I2 = normalize01(I)
%normalize01: normalize input vector, matrix, or graysale/color image pixel value to [0 1]
%   I2 = normalize01(I)
%
%Note: If all elements of I are equal, then all elements of I2 are 0.5;
%
% 2010.07.01 Xiang Huang Comments
%


Imax = double(max(I(:)));
Imin = double(min(I(:)));
if (Imax - Imin < 1e-6)
    %error('Input matrix has constant value, not able to normalize to [0, 1]');
    I2 = 0.5*ones(size(I));
    %I2 = I; %
else    
    I2 = (double(I) - Imin)/(Imax - Imin);  %normailze
end