function Sino2 = shiftForward(Sino, d, verbose)
%shiftForward:  Forward model to compute a sub-pixel shifted/drifted sinogram
% Sino:  NTheta * NTau, original sinogram
% d:     NTheta * 1, d(i) denotes the shift/drift amount of ith row of Sino
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
% Output:
%   Sino2: Shifted sinogram, each row n of the sino is (circularly) shifted to the right by d(n)
%         We expect user has zero-padded sino so that the circular shift will be linear-shift
[NTheta, NTau] = size(Sino);

Sino2 = zeros(NTheta, NTau, 'like', Sino);

for n = 1:NTheta
     % TODO, maybe vectorize it instead of for loop, but it's already fast
    %Sino2(n, :) = shiftForward1D(Sino(n,:), d(n)); % before 2020.02
    Sino2(n, :) = shiftForward1D(Sino(n,:), -d(n)); % after 2020.02, changed sign to match the direction of XTM_Tensor_XH()
    
end

if exist('verbose', 'var') && verbose
    if NTheta > 1
        subplot(121); imagesc(Sino); subplot(122); imagesc(Sino2);
    else
        plot(Sino, 'b'); hold on; plot(Sino2, 'r');
    end
end
end


%%
%
function f2 = shiftForward1D(f, d, verbose)
% (circular) shift along row to the right by amound d, which could be fractional
[~, Nc] = size(f);
getIdxFFT = @(N) ifftshift(-floor(N/2):ceil(N/2)-1); % notice the center of fft/ifft, need this fix
idxC =  getIdxFFT(Nc);
f2 = abs(ifft(fft(f).*exp(-1j*2*pi*d/Nc * idxC)));

if exist('verbose', 'var') && verbose
    plot(f, 'b'); hold on; plot(f2, 'r');
end
end


function f2 = shiftForward2D(f, dr, dc, verbose)
% (circular) shift of 2D along both row and column by amound d, which could be fractional
[Nr, Nc] = size(f);
getIdxFFT = @(N) ifftshift(-floor(N/2):ceil(N/2)-1); % notice the center of fft/ifft, need this fix
idxR = getIdxFFT(Nr);
idxC = getIdxFFT(Nc);
[idxC, idxR] = meshgrid(idxC, idxR);
f2 = abs(ifft2(fft2(f).*exp(-1j*2*pi*(dr*idxR/Nr + dc*idxC/Nc))));

if exist('verbose', 'var') && verbose
    subplot(121); imagesc(f); subplot(122); imagesc(f2);
end
end

function testMe(SAll, drift)
S1 = shiftForward(SAll(:,:,2), -drift);
S2 = shiftForward(SAll(:,:,1), drift);
difference(SAll(:,:,2), SAll(:,:,1))
difference(S2, SAll(:,:,2))
difference(S1, SAll(:,:,1))
figure(9999); multAxes(@imagesc, {SAll(:,:,1), SAll(:,:,2), S1, S2})
end


