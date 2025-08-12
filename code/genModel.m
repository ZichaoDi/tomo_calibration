function [S,L,WGT,LNormalizer] = genModel(sampleName,WGT,WSz,...
                                     NTheta,NTau, driftGT,LN)
% Input:
%   sampleName   - Specify a sample object, not case-sensitive.
%                  If WGT is nonempty, sampleName can be an empty string.
%                  [ 'Phantom' | 'Brain' | 'Golosio' | 'circle' |
%                   'checkboard' | 'fakeRod'| '' ]
%   WGT          - Groundtruth of sample object.
%                  If WGT needs to be generated, need to specify sampleName.
%   WSz          - Size of WGT.
%                  If WGT is nonempy, size(WGT) should be equal to WSz.
%   NTheta       - Number of rotations. Use odd number for display purposes.
%                  [ 30 | 45 | 60 | 90 | ...]
%   NTau         - Number of discrete beams. 
%                  Need to be large enough to cover object diagonal plus
%                  tolerance with maxDrift.
%   driftGT      - NTau by 1 vector of drifts for all beams.
%   noiseLevel   - Standard deviation of Gaussian noise added to sinogram.
%   LN           - Normalizer for L. 
%                  If nonzero, it is used directly;
%                  If zero, the function calculates a new LN.
%
% Output:
%   S            - NTheta by NTau sinogram.
%   L            - Forward Operator such that L*WGT(:) = S(:) if noisefree.
%   WGT          - Groundtruth of sample object if needs to be generated.
%   LNormalizer  - Normalizer for L if needs to be generated.


% maxDrift: maximum drift error; choose from [0.5, 1, 2, 4, 8] or [0.5, 1, 2, 3, 5] etc.
% NTau: % number of discrete beam, enough to cover object diagonal, 
% plus tolarence with maxDrift, also use + rem(NTau-Ny,2) to make NTau the 
% same odd/even as Ny just for perfection, so that for theta=0, we have 
% sum(WGT, 2)' and  S(1, (1:Ny)+(NTau-Ny)/2) are the same with a scale factor
% of sinagram of Phantom. As even NTheta will include theta of 90 degree where 
% sinagram will be very bright as Phantom sample has verical bright line on left and right boundary.

if isempty(WGT) && ~isempty(sampleName)
    WGT = createObject(sampleName, WSz); 
    WGT = WGT/max(WGT(:)); 
    assert(all(WGT(:)>=0), 'Groundtruth object should be nonnegative');
elseif isempty(WGT) && isempty(sampleName)
    error('Need to specify sampleName');
elseif size(WGT) ~= WSz
    error('WGT matrix should have size WSz');
end



L = XTM_Tensor_XH(WSz, NTheta, NTau, driftGT, WGT);
% want maximum row sum to be 1 rather than a too small number
if LN == 0
    LNormalizer = full(max(sum(L,2)));
else
    LNormalizer = LN;
end

L = L/LNormalizer;
S = reshape(L*WGT(:), NTheta, NTau);

end

