function [nnz_, nzMask] = nnz(x, relativeThreshold)
%[nnz_, nzMask] = compare.nnz(x, relativeThreshold); Compute the number/mask of non-zero elements whose absolute value > relativeThreshold.
% Input:
%    x:         nD numberic array
%   relativeThreshold:  relative relativeThreshold, positive, default is 10*eps(class(x))).
% Output:
%    nnz_:   the number of non-zero elements
%    nzMask: the binary mask with non-zero elements as true, same dimension as x (could be 1D/2D/3D/etc.)
% 2017.06.02:  Xiang Huang

if nargin < 2
   relativeThreshold = 10*eps(class(x)); 
else
    assert(relativeThreshold>0, 'Threshold=%.3e is smaller than 0, but it should be positive!');
end

nzMask = abs(x) > relativeThreshold*max(x(:));

nnz_ = sum(nzMask(:)); % equivalent but faster than nnz(abs(x(:)) > th)

end



