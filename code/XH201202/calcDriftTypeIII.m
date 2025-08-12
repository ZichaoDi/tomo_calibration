function drift = calcDriftTypeIII(Sm, S0_rec, maxDrift, interpMethod)
%compute the drift error from sinogram matching
% Input: 
%   S0_rec:  the Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
%   Sm: the Sinogram with beamline drifted, NTheta*NTau 2D matrix, assume it is interploation from S0_rec here. 
%   maxDrift: the maximum possible drift amount in the unit of beamlines, integer here (or pre-convert to integer: not tested)
%   interpMethod:    interplolation options
%               'linear'
%               'central_diff', 
% Output: 
%   drift: Type III drift. NTheta*NTau, assume small mostly [-0.5, 0.5], TODO: add check if it is too large
%          2020.07.06_still restrict drift to [-1, 1], todo change to [-0.5, 0.5]?
% See also
%   interpSinoTypeIII(); calcDrift();

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
% todo: for non-interger maxDrift, add further checking, add "if alphaCur + k < maxDrift"  before "if objCur <= objOpt"


%%
[NTheta, NTau] = size(S0_rec);
% 1st order difference: forward/central difference:
SDiffFwd = zeros(size(S0_rec),'like',S0_rec); n = 1:(NTau-1); SDiffFwd(:,n) = S0_rec(:,n+1) - S0_rec(:,n); % SDiffFwd = [S0_rec(:, 2:end) - S0_rec(:, 1:end-1), zeros(NTheta, 1)]; 


% drift = d+alpha, where d is the integer part, and 0 <= alpha <1 is the fractional part. 2020.02.22: limit d values to -1 or 0.
d = zeros(NTheta, NTau);
alpha = zeros(NTheta, NTau);

switch lower(interpMethod) 
    case lower('linear') 
        % Cases can't recover drift even if the Sm is generated from S0_rec with matching linear interpolation model
        % I. too smooth region, SDiffFwd(theta, tau2) = 0, this case is rare
        % II. two solutions, drift postive and negative both satisfy, this happens at a local miminum or maximum of sino
        % Case II has 50%*50% probablity wrong
        % to solve I & II, need to assume some proior on the drift, such from neighbors
        
         % loop version. Only keep for reference, same effect but slower than above vectorized version
        for theta = 1:NTheta
            for tau = 1:NTau
                %  theta=44, tau=35; drift(theta, tau)=0.0738, driftGT(theta,tau)=-0.0167
                if theta==44 && tau==35
                    % disp('debug');
                end
                objOpt = inf;
                for dCur = [-1, 0]
                    tau2 = tau + dCur; if tau2<1; tau2=1; end
                    %Sm(theta, tau) = S0_rec(theta, tau2) + SDiffFwd(theta, tau2)*alpha(theta, tau); %forward model for reference
                    %if SDiffFwd(theta, tau2) < eps % bug?
                    if abs(SDiffFwd(theta, tau2)) < eps
                        alphaCur = 0;                        
                    else
                        alphaCur = (Sm(theta, tau) - S0_rec(theta, tau2)) / SDiffFwd(theta, tau2);
                    end
                    alphaCur = min(max(alphaCur, 0), 1); % clip to [0, 1]. 
                    % if alphaCur == 1; alphaCur = alphaCur - 1e-3; end; % No need. When alphaCur = 1, next iteration will find alphaCur = 0 and overwrite it as alphaCur = 0. So we always have alphaCur in [0, 1)
                    objCur = abs(S0_rec(theta, tau2) + SDiffFwd(theta, tau2)*alphaCur - Sm(theta, tau));
                    %if objCur <= objOpt % bug, always pickup the later one if objCur == objOpt (== 0)
                    if (objCur < objOpt - eps) || (abs(objCur-objOpt)<=eps && abs(dCur+alphaCur) <= abs(d(theta, tau)+alpha(theta, tau)))
                        % condition after ||:  when objCur == objOpt, pickup the drift of the two has smaller absolute value
                        objOpt = objCur;
                        d(theta, tau) = dCur;
                        alpha(theta, tau) = alphaCur;
                    end
                end
            end            
        end
        
        % vectorized version. TODO
        %[Tau, Theta] = meshgrid(1:NTau, 1:NTheta);
        %Tau2 = Tau + d;  Tau2(Tau2<1) = 1;  Tau2(Tau2>NTau-1) = NTau-1; 
        % return the total drift
        drift = d + alpha;
    %case 'quadratic' 
    case 'central_diff'
        % test interpolation and reconstruction, assume drift is -1 to 1, and use pixels left and right for central difference
        % S0_rec is S, Sm is SInterp
        S = S0_rec; SInterp = Sm;
        SDiffCtr = zeros(size(S),'like',S); n = 2:(NTau-1); SDiffCtr(:,n) = (S(:,n+1) - S(:,n-1))/2; % SDiffCtr = [zeros(NTheta, 1), (S(:, 3:end) - S(:, 1:end-2))/2, zeros(NTheta, 1)];
        % SMeanCtr = zeros(size(S),'like',S); n = 2:(NTau-1); SMeanCtr(:,n) = (S(:,n+1) + S(:,n-1))/2; 
        SMeanCtr = S; 
        % SInterp = SDiffCtr.*drift + SMeanCtr; 
        drift = (SInterp - SMeanCtr) ./ SDiffCtr;
        drift(SDiffCtr == 0) = 0;
        drift = min(max(drift, -maxDrift), maxDrift); % clamp to [-maxDrift, maxDrift], NO need if SInterp = SDiffCtr.*drift + SMeanCtr is exact. 
    otherwise
        error('Unknown interpolation method: %s.', interpMethod);         
end




end

function testMe()
%% 2020.03_tested
% 2020.06.03_tested again for 50x50 Phantom with sinogram 45x74, drift std_=0.5
S0_rec = SAll(:,:,1);
mask = abs(S0_rec) > eps;

%% I. Test use drifted sinogram generated from interpSinoTypeIII()
S2_drifted_via_interp = interpSinoTypeIII(S0_rec, driftGT, 'central_diff');  % 'linear' (default), and 'central_diff'
drift = calcDriftTypeIII(S2_drifted_via_interp, S0_rec, 1,  'central_diff');  % 'linear' (default), and 'central_diff'
figure; showDriftDifference(drift, driftGT, mask);

%% II. Test use drifted sinogram generated from measurement (with exact forward model)
Sm = SAll(:,:,2);
drift = calcDriftTypeIII(Sm, S0_rec, 1, 'central_diff');  % 'linear' (default), and 'central_diff'
figure; showDriftDifference(drift, driftGT, mask);

%% III. Test use drifted sinogram generated from exact forward model with super-resolution drift grid
% TODO: optimize to reduce computation by half, only shift positive, as shift negative could be computed by shift left neighbor
% positive
N_res = 10; %2, 5, 10
%d_grid = (0:N_res-1)/N_res;
%drift_resAll = bsxfun(@plus, reshape(d_grid, N_res, 1, 1), zeros(1, NTheta, NTau));
d_grid = (-(N_res-1):N_res-1)/N_res;
N_grid = numel(d_grid);
drift_resAll = bsxfun(@plus, zeros(NTheta, NTau, 1), reshape(d_grid, 1, 1, N_grid));

L_resAll = cell(N_grid, 1);
for n = 1:N_grid
    L_resAll{n} = XTM_Tensor_XH(WSz, NTheta, NTau, drift_resAll(:, :, n));        
end
for n = 1:N_grid
    L_resAll{n} = L_resAll{n} / LNormalizer; % LNormalizer
end

S_resAll = zeros(NTheta, NTau, N_grid);
for n = 1:N_grid
    S_resAll(:, :, n) =  reshape(L_resAll{n}*WGT(:), NTheta, NTau);
end
% hack, flip to get -1, uses more memory to simplify coding for now

Sm = SAll(:,:,2);

%%
% compute drift by comparing S2_drfited_via_exactforwardmodel and S_resAll
err = abs(bsxfun(@minus, Sm, S_resAll));
%drift = calcDriftTypeIII(S2_drfited_via_exactforwardmodel, S0_rec, 1,  'central_diff');  % 'linear' (default), and 'central_diff'
drift = zeros(NTheta, NTau);
[~, idx] = min(err,[],3); 
for n = 1:numel(idx)
    drift(n) = d_grid(idx(n));
end
drift(Sm < eps) = 0;

figure; showDriftDifference(drift, driftGT, mask);

end
