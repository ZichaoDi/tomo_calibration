function SInterp = interpSino(S, drift, interpMethod)
%INTERPSINO polynomial interpolation of sinogram
% Input:   
%       S:      2D original sinogram of size NTau*NTheta
%       drift:  the drift amount in object pixel units for each beamline of size NTau*1. We assume that each beam has same drift for different NTheta
%       interpMethod:    interplolation options
%               'linear'
%               'quadratic', 
% Output:   
%       SInterp: 2D interpolated sinogram of size NTau*NTheta
% See also:
%       calcDrift(); genDriftMatrix();

if nargin < 1; testMe(); end

if ~exist('interpMethod', 'var') || isempty(interpMethod)
    interpMethod = 'linear';
end


[NTheta, NTau] = size(S);
SInterp = zeros(size(S), 'like', S);
% 1st order difference: forward/central difference
SDiffFwd = zeros(size(S),'like',S); n = 1:(NTau-1); SDiffFwd(:,n) = S(:,n+1) - S(:,n); % SDiffFwd = [S(:, 2:end) - S(:, 1:end-1), zeros(NTheta, 1)]; 
SDiffCtr = zeros(size(S),'like',S); n = 2:(NTau-1); SDiffCtr(:,n) = (S(:,n+1) - S(:,n-1))/2; %SDiffCtr = [zeros(NTheta, 1), (S(:, 3:end) - S(:, 1:end-2))/2, zeros(NTheta, 1)];
%  2nd order difference
SDiff2nd = zeros(size(S),'like',S); n = 2:(NTau-1); SDiff2nd(:,n) = (S(:,n+1) + S(:,n-1) - 2*S(:,n)); % SDiff2nd = [zeros(NTheta, 1), S(:, 3:end) + S(:, 1:end-2) - 2*S(:,2:end-1), zeros(NTheta, 1)]; 


% compute integer and fractional part, which is used for linear/quadratic opitions
d = floor(drift);  alpha = drift - d;  % 0 <= alpha < 1

switch lower(interpMethod)
    case {lower('linearRaw'), 'linear_raw'}
        % we assume -1< = drift <=1, we use beam n-1 and n+1 (of S) to recover beam n (of SInterp).
        % This is only used as template and comparison, not as good as linear interpolation that restrict to  0<= alpha <1, 
        % SInterp(:,n) = S(:,n) + drift(n)*(S(:,n+1) - S(:,n-1))/2; % 1st order Taylor expansion with central difference
        SInterp = S + bsxfun(@times, SDiffCtr, drift');       
    case {lower('quadraticRaw'),'quadratic_raw'}
        % quaractic interpolation:  assume -1< = drift <=1, we use beam n-1, n, and n+1 (of S) to recover beam n (of SInterp).
        % This is only used as template and comparison, not as good as quadratic interpolation that restrict to  -0.5<= alpha <0.5, 
        % SInterp(:,n) = S(:,n) + drift(n)*(S(:,n+1) - S(:,n-1))/2 + 0.5*drift(n)^2*(S(:,n+1) + S(:,n-1) - 2*S(:,n))  ; %2nd order Taylor expansion      
        SInterp = S + bsxfun(@times, SDiffCtr, drift') +  bsxfun(@times, SDiff2nd, 0.5*drift'.^2); 
    case 'linear'
        % To recover beam i of SInterp, we use two cloest beams j and j+1 (of S) to recover beam i (of SInterp), where j = beams(i).        
        % SInterp(:, i) =  S(:, j) + alpha(i)*(S(:,j+1) - S(:,j)); % 0 <= alpha(i) < 1 is the drifted distance to beam j (in S) from beam i (in SInterp)        
        beams = (1:NTau)' + d;  beams(beams<1) = 1;  beams(beams>NTau) = NTau;  % beams(i), beams(i)+1 are the two closest beamline in S for nth beamline in SInterp
        SInterp = S(:, beams) + bsxfun(@times, SDiffFwd(:, beams), alpha');         
    case 'quadratic'
        % quaractic interpolation similar to quadraticRaw, but restrict -0.5 <= alpha <= 0.5, where drift = d + alpha and d is interger       
        % To recover beam i of SInterp, we use beam 
        idx = alpha >= 0.5;
        d(idx) = d(idx) + 1; alpha(idx) = alpha(idx) - 1;  % -0.5 <= alpha(n) < 0.5 is the drifted distance to beams(n) (of S) for beam n (of SInterp)
        beams = (1:NTau)' + d;  beams(beams<1) = 1;  beams(beams>NTau) = NTau;  
        SInterp = S(:, beams) + bsxfun(@times, SDiffCtr(:, beams), alpha') +  bsxfun(@times, SDiff2nd(:,beams), 0.5*alpha'.^2);         
    case {lower('linearOld'), 'linear_old'}
        % old implementation with same result as 'linear', but more complicated and slower, and boundary condition bad, just leave for comparison
        for i = 1:NTau
            j = i+d(i);
            %SInterp(:, i) = (1-alpha(i)) * S(:, j)  + alpha(i)*S(:, j+1);
            if (j >= 1) && (j <= NTau)
                SInterp(:, i) = (1-alpha(i)) * S(:, j);
                if (alpha(i) > eps)
                    if ((j+1) <= NTau)  % j >= 1 already, no need to check j+1 >= 1
                        SInterp(:, i) = SInterp(:, i) + alpha(i)*S(:, j+1);
                    else
                        error('Out of boundary! You should not reach here: %d, alpha(i)=%.2f.\n', d(i), alpha(i));
                    end
                end
                %NOTE: if alpha(i) < eps, then S(:, j+1) is not used, and we no need to check whether j+1 is within [1, NTau].
            else
                error('Out of boundary! You should not reach here: %d, alpha(i)=%.2f.\n', d(i), alpha(i));
            end
        end        
    otherwise
        error('Unknown interpolation method: %s.', interpMethod);        
end



end


function testMe()
%% Test the quality of different interpolation method
load('ret_20181121_phantom.mat'); %SAll, X, driftAll, driftGT, maxDrift
interpMethods = {'linear', 'quadratic', 'linearRaw', 'quadraticRaw', 'linearOld'}; % [34.82, 34.99, 30.41, 34.74, 34.82]dB
for n = 1:numel(interpMethods)
    S2All{n} = interpSino(SAll(:,:,1), driftGT, interpMethods{n});
end
for n = 1:numel(interpMethods)    
    fprintf('Interpolation method: %s, psnr = %.2fdB.\n', interpMethods{n}, difference(S2All{n},  SAll(:,:, 2)));   
end
%%
figure(1); clf; multAxes(@imagesc, {S2All{:}, SAll(:,:,2)}); multAxes(@title, {interpMethods{:}, 'original'}); multAxes(@colorbar); linkAxesXYZLimColorView(); 
%%
%% Reconstruction with solver from XH, with L1/TV regularizer.
%{
% Need 100/500/1000+ iteration to get good/very good/excellent result with small regularizer.
% choose small maxSVDofA to make sure initial step size is not too small. 1.8e-6 and 1e-6 could make big difference for n=2 case
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
regWt = 1e-6; % 1e-6 to 1e-8 both good, b
maxIterA = 500; % 100 is not enough
maxSVDofA = 1e-6; %svds(L, 1)*1e-4; % multiply by 1e-4 to make sure it is small enough so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA};

WInterp = zeros(WSz(1), WSz(2), numel(S2All), 'like', WGT);
for n = 1:numel(S2All)
    WInterp(:,:,n) = solveTwist(S2All{n}, LAll{n}, paraTwist{:});
end
figNo = 500;
for n = 1:numel(S2All)    
    figure(figNo+n); W = WInterp(:,:,n);  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', interpMethods{n}, difference(W, WGT), regType, regWt, maxIterA));
end

figNo = 500;
for n = 1:numel(S2All)    
    figure(figNo+n); W = WInterp(:,:,n);  imagesc(W-WGT); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', interpMethods{n}, difference(W, WGT), regType, regWt, maxIterA));
end
tilefigs;
%}

end
