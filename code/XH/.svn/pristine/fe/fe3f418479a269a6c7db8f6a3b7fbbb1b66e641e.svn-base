function d = shiftBackward(Sino, Sino2, usfac, verbose)
%shiftBackward:  Back model to compute a sub-pixel shift between two sinograms, from Sino to Sino2
%   Sino:  NTheta * NTau, original sinogram
%   Sino2: Shifted sinogram, each row n of the sino is (circularly) shifted to the right by d(n)
%   usfac: Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 10)
% Output:
%   d:     NTheta * 1, d(i) denotes the shift/drift amount of ith row of Sino
%         We expect user has zero-padded sino so that the circular shift will be linear-shift
[NTheta, ~] = size(Sino);
if ~exist('usfac', 'var') || isempty(usfac)
    usfac = 10;
end
d = zeros(NTheta, 1);

for n = 1:NTheta    
    [output, ~] = dftregistration(fft2(Sino(n,:)),fft2(Sino2(n,:)), usfac);
    %d(n) = -output(4); % before 2020.02, 
    d(n) = output(4); % after 2020.02, changed sign to match the direction of XTM_Tensor_XH()
end

if exist('verbose', 'var') && verbose
    plot(d);
end
end


%%
%
% function g = shiftBackward1D(f, f2, verbose)
% end
% 
% function g = shiftBackward2D(f, f2, verbose)
% end
% 

function testMe(SAll, driftGT)
drift = shiftBackward(SAll(:,:,1), SAll(:,:,2));
figure; plot(driftGT, 'b'); hold on; plot(drift, 'r'); legend('drift-GT', 'drift-Rec');
difference(drift, driftGT)
drift_rev = shiftBackward(SAll(:,:,2), SAll(:,:,1));
difference(-drift_rev, drift)
end