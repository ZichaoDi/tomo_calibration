function drift = calcDrift(S2, S, maxDrift, interpMethod)
%compute the drift error from sinogram matching
% Input: 
%   S2:  the measured Sinogram with beamline drifted, NTheta*NTau 2D matrix, assume it is interploation from S here. 
%   S: the estimated Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
%   maxDrift: the maximum possible drift amount in the unit of beamlines, integer here (or pre-convert to integer: not tested)
%   interpMethod:    interplolation options
%               'linear'
%               'quadratic', 
% Output: 
%   drift: the drift amount for each beamline (assume invariant with rotation angle),  NTau*1 vector
%               TODO?: change to 1*NTau for consistency of other cases: NTheta*1 and NTheta*NTau drifts?
% See also
%   interpSino(); genDriftMatrix();
%  2019.07.10: revise linear model for same notation of a, b
%  2019.07.11-08.15: quadratic model done


if nargin < 1; testMe(); end


if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

if ~exist('interpMethod', 'var') || isempty(interpMethod)
    interpMethod = 'linear';
end

%assert(maxDrift==floor(maxDrift), 'maxDrift need to be integer!, but it is %.5f. TODO: correct this', maxDrift)
maxDriftCeil = ceil(maxDrift);
if maxDrift ~=  maxDriftCeil
    warning('non-integer maxDrift %.5f not fully tested', maxDrift)
end
% todo: for non-interger maxDrift, add further checking, add "if alphaCur + dCur < maxDrift"  before "if objCur <= objOpt"


%%
[NTheta, NTau] = size(S2);
% 1st order difference: forward/central difference:
SDiffFwd = zeros(size(S),'like',S); n = 1:(NTau-1); SDiffFwd(:,n) = S(:,n+1) - S(:,n); % SDiffFwd = [S(:, 2:end) - S(:, 1:end-1), zeros(NTheta, 1)]; 
SDiffCtr = zeros(size(S),'like',S); n = 2:(NTau-1); SDiffCtr(:,n) = (S(:,n+1) - S(:,n-1))/2; % SDiffCtr = [zeros(NTheta, 1), (S(:, 3:end) - S(:, 1:end-2))/2, zeros(NTheta, 1)];
%  2nd order difference
SDiff2nd = zeros(size(S),'like',S); n = 2:(NTau-1); SDiff2nd(:,n) = (S(:,n+1) + S(:,n-1) - 2*S(:,n)); % SDiff2nd = [zeros(NTheta, 1), S(:, 3:end) + S(:, 1:end-2) - 2*S(:,2:end-1), zeros(NTheta, 1)]; Same as SDiff2nd(:,n) = SDiffFwd(:,n) - SDiffFwd(:,n-1);



% drift = d+alpha, where d is the integer part, and 0 <= alpha <1 is the fractional part
d = zeros(NTau, 1);
alpha = zeros(NTau, 1);

switch lower(interpMethod)
    case lower('linearOld')  % old a, b notations before 2019.07.11, only good for interger maxDrift
        for i = maxDriftCeil+1 : NTau-maxDriftCeil % compute d(i), note zero boundaries d(i)=0 for i < maxDriftCeil+1 or i > NTau-maxDriftCeil
            if norm(S2(:,i), 2) < NTheta*eps && norm(S(:,i), 2) < NTheta*eps; continue;  end  % empty beam matched to empty beam, skip & set drift to 0
            objOpt = inf;
            for dCur = -maxDriftCeil : maxDriftCeil-1  % interpolate between beam i0 and i0+1, where i0 = i+dCur;
                i0 = i+dCur; 
                % find the optimal alpha(i) that minimize ||b - a*alpha ||^2
                a = SDiffFwd(:, i0);
                b = S2(:,i) - S(:,i0);
                aTa = a'*a;
                aTb = a'*b;
                bTb = b'*b;
                if aTa < NTheta*eps
                    alphaCur = 0;
                else
                    alphaCur = aTb / aTa;
                end
                alphaCur = min(max(alphaCur, 0), 1); % clip to [0, 1]. 
                % if alphaCur == 1; alphaCur = alphaCur - 1e-3; end; % No need. When alphaCur = 1, next iteration will find alphaCur = 0 and overwrite it as alphaCur = 0. So we always have alphaCur in [0, 1)
                objCur = aTa*alphaCur^2 - 2*aTb*alphaCur + bTb;
                if objCur <= objOpt
                    objOpt = objCur;
                    d(i) = dCur;
                    alpha(i) = alphaCur;
                end
            end
        end
    case lower('linear') % new a, b notations after 2019.07.12, only good for interger maxDrift
        for i = maxDriftCeil+1 : NTau-maxDriftCeil % compute d(i), note zero boundaries d(i)=0 for i < maxDriftCeil+1 or i > NTau-maxDriftCeil
            if norm(S2(:,i), 2) < NTheta*eps && norm(S(:,i), 2) < NTheta*eps; continue;  end  % empty beam matched to empty beam, skip & set drift to 0
            objOpt = inf;
            for dCur = -maxDriftCeil : maxDriftCeil-1  % interpolate between beam i0 and i0+1, where i0 = i+dCur;
                i0 = i+dCur; 
                % find the optimal alpha(i) that minimize ||a + b*alpha ||^2
                a = S(:,i0)-S2(:,i);
                b = SDiffFwd(:, i0);
                bTb = b'*b;
                bTa = b'*a;
                aTa = a'*a;
                if bTb < NTheta*eps
                    alphaCur = 0;
                else
                    alphaCur = - bTa / bTb;
                end
                alphaCur = min(max(alphaCur, 0), 1); % clip to [0, 1]. 
                % if alphaCur == 1; alphaCur = alphaCur - 1e-3; end; % No need. When alphaCur = 1, next iteration will find alphaCur = 0 and overwrite it as alphaCur = 0. So we always have alphaCur in [0, 1)
                objCur = bTb*alphaCur^2 + 2*bTa*alphaCur + aTa;
                if objCur <= objOpt
                    objOpt = objCur;
                    d(i) = dCur;
                    alpha(i) = alphaCur;
                end
            end
        end
    case lower('quadraticTest')
        % don't loop through all dCur, but rather compute and store result for all posible dCur, and pickup the best        
        for i = maxDriftCeil+1 : NTau-maxDriftCeil % compute d(i), note zero boundaries d(i)=0 for i < maxDriftCeil+1 or i > NTau-maxDriftCeil
            if norm(S2(:,i), 2) < NTheta*eps && norm(S(:,i), 2) < NTheta*eps; continue;  end  % empty beam matched to empty beam, skip & set drift to 0
            objOpt = inf;
            for dCur = -maxDriftCeil : maxDriftCeil  % interpolate between beam i0-1, i0, and i0+1, where i0 = i+dCur. Note the difference from 'linear' where we loop dCur from -maxDriftCeil to maxDriftCeil, not to maxDriftCeil-1 
                % find the optimal alpha(i) that minimize f(alpha) = ||a + b*alpha + c*alpha^2 ||^2
                i0 = i+dCur;
                a = S(:,i0) - S2(:,i);
                b = SDiffCtr(:, i0);
                c = SDiff2nd(:, i0);
                
                % f(alpha) = pObj.*[alpha^4, alpha^3, alpha^2, alpha^1, 1]                
                pObj = [(c'*c)/4, c'*b, (b'*b + c'*a), 2*(b'*a), a'*a];
                % f'(alpha) = pGrad.*[alpha^3, alpha^2, alpha^1, 1]
                pGrad = pObj(1:end-1) .* [4, 3, 2, 1];
                % find the exterme points that makes pGrad = 0
                alphaCurAll = roots(pGrad);
                
                if ~isempty(alphaCurAll)    % alphaCurAll could be empty for special cases such as pGrad = [0, 0, 0, 0]
                    % only keep real alpha that fall in (-0.5, 0.5)
                    alphaCurAll = alphaCurAll(isreal(alphaCurAll));
                    alphaCurAll = alphaCurAll(alphaCurAll > -1/2 & alphaCurAll < 1/2);                
                    % add boundary points to candidate points
                    alphaCurAll = [alphaCurAll, -0.5, 0.5]; % todo: maybe only add -0.5?, as 0.5 will be equivalent to -0.5 in next iteration
                    % make sure those points satisfy dCur+alpha within [-maxDrift, maxDrift]
                    alphaCurAll = alphaCurAll(dCur+alphaCurAll > -maxDrift & dCur+alphaCurAll < maxDrift);
                end
                
                % pickup the best points, if the candidates set is not empty
                if ~isempty(alphaCurAll)
                    objCurAll = polyval(pObj, alphaCurAll);
                    [objCur, idx] = min(objCurAll);
                    alphaCur = alphaCurAll(idx);

                    % compare with previous best computed drift(i)
                    if objCur <= objOpt
                        objOpt = objCur;
                        d(i) = dCur;
                        alpha(i) = alphaCur;
                    end 
                end
            end
            assert(objOpt ~= inf, 'objOpt inf, means not finding any alpha?');
        end            
    case 'quadratic'
        % don't loop through all dCur, but rather compute and store result for all posible dCur, and pickup the best        
        % todo: check why so many outputs of alpha = 0.5 or -0.5
        for i = maxDriftCeil+1 : NTau-maxDriftCeil % compute d(i), note zero boundaries d(i)=0 for i < maxDriftCeil+1 or i > NTau-maxDriftCeil
            if i == 118
                fprintf('debug\n');
            end
            if norm(S2(:,i), 2) < NTheta*eps && norm(S(:,i), 2) < NTheta*eps; continue;  end  % empty beam matched to empty beam, skip & set drift to 0
            objOpt = inf;
            for dCur = -maxDriftCeil : maxDriftCeil  % interpolate between beam i0-1, i0, and i0+1, where i0 = i+dCur. Note the difference from 'linear' where we loop dCur from -maxDriftCeil to maxDriftCeil, not to maxDriftCeil-1 
                % find the optimal alpha(i) that minimize f(alpha) = ||a + b*alpha + c*alpha^2 ||^2
                i0 = i+dCur;
                a = S(:,i0) - S2(:,i);
                b = SDiffCtr(:, i0);
                c = SDiff2nd(:, i0);
                
                % f(alpha) = pObj.*[alpha^4, alpha^3, alpha^2, alpha^1, 1]                
                pObj = [(c'*c)/4, c'*b, (b'*b + c'*a), 2*(b'*a), a'*a];
                % f'(alpha) = pGrad.*[alpha^3, alpha^2, alpha^1, 1]
                pGrad = pObj(1:end-1) .* [4, 3, 2, 1];
                % find the exterme points that makes pGrad = 0
                alphaCurAll = roots(pGrad)'; % col vector to vector
                
                if ~isempty(alphaCurAll)    % alphaCurAll = [] for special cases such as pGrad = [0, 0, 0, 0]; but isreal([])=1
                    % only keep real alpha
                    alphaCurAll = alphaCurAll(abs(imag(alphaCurAll)) < eps); % 
                    alphaCurAll = real(alphaCurAll);
                    % only keep alpha in (-0.5, 0.5)
                    alphaCurAll = alphaCurAll(alphaCurAll > -1/2 & alphaCurAll < 1/2);                
                end
                % add boundary points to candidate points
                alphaCurAll = [alphaCurAll, -0.5, 0.5]; % todo: maybe only add -0.5?, as 0.5 will be equivalent to -0.5 in next iteration
                % make sure those points satisfy dCur+alpha within [-maxDrift, maxDrift]
                alphaCurAll = alphaCurAll(dCur+alphaCurAll > -maxDrift & dCur+alphaCurAll < maxDrift);
                
                % pickup the best points, if the candidates set is not empty
                if ~isempty(alphaCurAll)
                    objCurAll = polyval(pObj, alphaCurAll);
                    [objCur, idx] = min(objCurAll);
                    alphaCur = alphaCurAll(idx);

                    % compare with previous best computed drift(i)
                    if objCur <= objOpt
                        objOpt = objCur;
                        d(i) = dCur;
                        alpha(i) = alphaCur;
                    end 
                end
            end
            assert(objOpt ~= inf, 'objOpt inf, means not finding any alpha?');
        end      
    otherwise
        error('Unknown interpolation method: %s.', interpMethod);         
end

% return the total drift
drift = d + alpha;


end

function testMe()
%%
%clear; load('ret_20181121_phantom.mat'); % SAll, X, driftAll, driftGT, maxDrift % [20.09, 22.10 ]dB for {'quadratic', 'linear'}, quaratic bad due to boundry beams? to check
%clear; load('20190711_phantom.mat'); % SAll, driftGT, maxDrift, WGT, LAll, guassianSTD % [17.12, 19.27 ]dB for  {'quadratic', 'linear'}
clear; load('20190711_brain.mat'); % SAll, driftGT, maxDrift, WGT, LAll, guassianSTD % [26.23, 26.06 ]dB for {'quadratic','linear'}, quaratic bad due to boundry beams? to check
%maxDrift = maxDrift + 1e-8; % test non-interger maxDrift
drift = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, 'linear');

interpMethods = {'quadratic', 'linear'}; % [17.12, 19.27, 19.27]dB
NMethods = numel(interpMethods);
driftAll = zeros(numel(driftGT), NMethods);
for n = 1:NMethods
    driftAll(:,n) = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, interpMethods{n});
end

% Clean the left/right border, set to 0
drift = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, 'linear');
driftGT2 = driftGT;
iStart = find(cumsum(drift), 1, 'first');
iEnd = numel(driftGT) + 1 - find(cumsum(drift(end:-1:1)), 1, 'first');
driftGT2([1:iStart-1, iEnd+1:end]) = 0;
%
psnrAll = zeros(1, NMethods);
for n = 1:NMethods
    psnrAll(n)=difference(driftAll(:,n),  driftGT2);
    fprintf('Interpolation method: %s, psnr = %.2fdB.\n', interpMethods{n}, psnrAll(n));   
    figure(10+n); clf; plot([driftGT2, driftAll(:,n)], '-x'); title(sprintf('psnr=%.2fdB', psnrAll(n))); legend('GT', interpMethods{n});
end
fprintf('-------'); fprintf('%.2f, ', psnrAll); fprintf('-------\n'); 
%%
fprintf('verify that we really find the alpha that minimize the drift. Sinogram Interpolation from driftCalc should match actual drifted sinogram better than from drifGT. \n');
for n = 1:numel(interpMethods)
    S2All_driftGT{n} = interpSino(SAll(:,:,1), driftGT2, interpMethods{n});
    S2All_driftCalc{n} = interpSino(SAll(:,:,1), driftAll(:,n), interpMethods{n});
end
for n = 1:numel(interpMethods)    
    fprintf('Interpolation method from driftGT: %s, psnr = %.2fdB.\n', interpMethods{n}, difference(S2All_driftGT{n},  SAll(:,:, 2)));   
    fprintf('Interpolation method from driftCalc: %s, psnr = %.2fdB.\n', interpMethods{n}, difference(S2All_driftCalc{n},  SAll(:,:, 2)));   
end
fprintf('-------\n'); 



end
