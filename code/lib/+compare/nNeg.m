function [nNeg_, negMask] = nNeg(x, threshold)
%[nNeg_, negMask] = compare.nNeg(x, threshold); Compute the number/mask of negative elements whose value < -threshold. 
% Input:
%    x:         nD numberic array
%   threshold:  positive threshold, default is 10*eps(class(x))).
% Output:
%    nNeg_:   the number of negative elements.
%    negMask: the binary mask with negative elements as true, same dimension as x (could be 1D/2D/3D/etc.)
% 2017.06.01:  Xiang Huang


if nargin < 2
   threshold = 10*eps(class(x)); 
else
    assert(threshold>0, 'Threshold=%.3e is smaller than 0, but it should be positive!');
end

negMask = x < -threshold;

nNeg_ = sum(negMask(:)); % equivalent but faster than nnz(abs(x(:)) > th)

end

