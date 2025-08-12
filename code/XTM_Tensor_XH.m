function [L, DisR_Simulated] = XTM_Tensor_XH(WSz, NTheta, NTau, drift, W, n_delta, ind_cr)
% todo add W
%%%Simulate XRT of a given object with given detector and beam
% Inputs:
%   WSz:          The sample object size: 2*1 array defines [Ny, Nx]
%   NTheta, NTau: The numberical resolution of Sinogram. We may also use a single variable SSz = [NTheta, NTau] or [NTau, NTheta]
%   drift:        Optional drift amount of each beam. NTau*1 array.  Assume the drift is the same for all the angles for the same beam
%   W:            The optional sample object Ny*Nx image. We don't need W to generate the forward model L, but need it to
%                 generate DisR_Simulated and show the plots
%   n_delta:      controls center of rotation errors, default is 0 where there is no error.
%   ind_cr:       controls center of rotation errors, default set to 1 where there is no error
% Output:
%   L:             the forward model from sample object to tomography
%   DisR_Simulated: 
%
% 2018.09.04:   XH replaced nTau+1 --> NTau (strange and easily confused with the +1), NTheta ---> NTheta
% 2018.11.30:   XH moved generate object function out, As we may scan multiple times with different drift for the same object.
% 2018.12.07:   XH made the code much faster! Speed up of x3.4 times for middle scale problem and maybe more for larger problems.
%               For example, old code takes 25.4 seconds, new code takes 7.5 seconds for WSz = [100, 100]; NTheta = 90; NTau = 144; 
%               Old code called almost NTheta*NTau times of:  L(sub2ind([NTheta,NTau],n,i), currentInd) = Lvec;% /sum(Lvec);
%               Each time, matlab may need to keep reallocate memory for L matrix because the number of nonzeros are changing.
%               New cold: just assemable sparse matrix once in the end:  L = sparse(i,j,s,m,n);


if ~exist('drift', 'var') || isempty(drift)
    drift = zeros(NTau, 1);
end

if ~exist('W', 'var') || isempty(W)
    isWDefined = 0;
else
    isWDefined = 1;
    assert(all(size(W)==WSz), 'The given size and the actual size of sample object does not match!');
end

if ~exist('n_delta', 'var') || isempty(n_delta)
    n_delta = 0;
end
if(~exist('ind_cr','var'))
    ind_cr=1;
end

% user options to display scan results related to a given sample W
plotDisBeam =  0;  % display beam movement 
plotRotation = 0;  % display and debug for COR drift?
if ~isWDefined
    % force no display when W is not given
    plotDisBeam = 0;
    plotRotation = 0; 
end
if(plotDisBeam || plotRotation)
    figure(9999); clf;
end


more off;
% physical unit of object
% Tol = 1e-2; % ??
Tol = 1e-3;
omega = [-1,  1,  -WSz(1)/WSz(2),  WSz(1)/WSz(2)].*Tol; %omega = [-1     1    -1     1].*Tol; % XH: [xMin, xMax yMin yMax] in physical unit
dxdy = [(omega(2)-omega(1))/WSz(2) (omega(4)-omega(3))/WSz(1)]; % pixel size for [x,y] in physical unit? dx/dy
assert( abs(dxdy(1)/dxdy(2)-1) < eps, 'Currently we need to assume object has square pixel size dx = dy!');

% set delta_d0, delta0: the center of rotaion (COR) drift: error, default no error, so set delta_d0, delta0 to all zero
[delta_d0, delta0] = getCORDrift(n_delta, ind_cr, WSz, NTheta, dxdy);

% Setup the beam source and Detector
[Theta, SourceKnot0, DetKnot0] = Define_Detector_Beam_Gaussian_XH(omega, dxdy, delta_d0, NTheta, NTau);

%%==============================================================
NonEmptyBeam=[];

% L = sparse(NTheta*NTau, prod(WSz));


xbox=[omega(1) omega(1) omega(2) omega(2) omega(1)];  % XH: [xMin xMin xMax xMax xMin], xbox & ybox defined clockwise trasveral order from bottom-left corner
ybox=[omega(3) omega(4) omega(4) omega(3) omega(3)];  % XH: [yMin yMax yMax yMin yMin]


SourceKnot0(:,2) = SourceKnot0(:,2) + drift*dxdy(1);
DetKnot0(:,2) = DetKnot0(:,2) + drift*dxdy(1);

% Discretize the whole image space and get the boundary coordinates of WSz(2)*WSz(1) pixels
x = linspace(omega(1), omega(2), WSz(2)+1);  % linspace(xMin, xMax, Nx+1)
y = linspace(omega(3), omega(4), WSz(1)+1);

ii = []; jj = []; vv = [];  % XH: L = sparse(ii,jj,vv, NTheta*NTau, prod(WSz)) generates a sparse matrix L of size NTheta*NTau by prod(WSz) from the triplets ii, jj, and vv such that L(ii(k),jj(k)) = vv(k).

% XH: found a bug for 90 degree
for n=1:NTheta
    if (n==30 || n==31)
        1; % XH: found a bug for 90 degree when NTheta=60, n=31
    end
    theta = Theta(n); % current angle in radians
    thetaD = theta*180/pi; % current angle in degrees
    if(mod(n,10)==0)
        fprintf(1,'====== Angle Number  %d of %d: %d%c\n', n, NTheta, thetaD, char(176));
    end
    % rotate clockwise by theta (i.e. counter-clockwise by -theta), plus error delta0 (drift of center of rotation) which is 0 for XH
    TransMatrix = full([cos(theta), sin(theta), delta0(n,1)-cos(theta)*delta0(n,1)-sin(theta)*delta0(n,2); ...
                        -sin(theta), cos(theta), delta0(n,2)+sin(theta)*delta0(n,1)-cos(theta)*delta0(n,2);...
                        0 0 1]); 
    DetKnot = TransMatrix*[DetKnot0'; ones(1,size(DetKnot0,1))];
    DetKnot = DetKnot(1:2,:)';
    SourceKnot = TransMatrix*[SourceKnot0'; ones(1,size(DetKnot0,1))];
    SourceKnot = SourceKnot(1:2,:)';
    %% =========================================
    if(plotDisBeam)
        hold on;                
        imagesc([omega(1),omega(2)]+[0.5, -0.5]*dxdy(1), [omega(3),omega(4)]+[0.5, -0.5]*dxdy(2), W);
        title(sprintf('%.1f%c', thetaD, char(176)));%title(sprintf([num2str(int16(thetaD)),'%c'], char(176)));        
        plot(SourceKnot(:,1),SourceKnot(:,2),'r.-',DetKnot(:,1),DetKnot(:,2),'k.-');
        axis equal;
        axis([-1, 1, -1, 1]/2 * sqrt(2)*(max(SourceKnot0(:,2))-min(SourceKnot0(:,2))));
        %pause;
        drawnow;
    end
    if(plotRotation)
        output(:,:,:,n) = rotateAround(padarray(W,[ceil(WSz(1)/2),ceil(WSz(1)/2)]),delta0(n,1)/dxdy(1)+WSz(1),delta0(n,2)/dxdy(1)+WSz(1), thetaD);
        % subplot(5,ceil(NTheta/5),n)
        imagesc(sum(output(:,:,:,n),3)); axis xy image; hold on; plot(delta0(n,1)/dxdy(1)+WSz(1),delta0(n,2)/dxdy(1)+WSz(1),'r*');
        hold on;
        set(gca,'xtick',[],'ytick',[]);
        drawnow;
    end
    for i=1:NTau
        %============================= Plot Grid and Current Light Beam =================
        [index, Lvec, linearInd] = IntersectionSet_XH(SourceKnot(i,:), DetKnot(i,:), xbox, ybox, theta, ...
            x, y, omega, WSz, dxdy, Tol);
        %%%%%%%%================================================================
        if(~isempty(index) && norm(Lvec)>0)
            NonEmptyBeam = [NonEmptyBeam, sub2ind([NTheta,NTau],n,i)];
            currentInd = sub2ind(WSz,index(:,2),index(:,1));
            %L(sub2ind([NTheta,NTau],n,i), currentInd) = Lvec;% /sum(Lvec); % Old code! bad
            %XH: todo: fix the array size ii in the beginning. It should be less than NTheta*NTau*WSz? but this speedup is minor
            %as ii, jj, vv are just 1D array and matlab or any modern complier may just use dynamic array allocation technology
            ii = [ii; repmat(sub2ind([NTheta,NTau],n,i), [numel(Lvec), 1])]; 
            jj = [jj; currentInd]; 
            vv = [vv; Lvec]; 
        end
    end

end

L = sparse(ii,jj,vv, NTheta*NTau, prod(WSz)); % generates a sparse matrix L of size NTheta*NTau by prod(WSz) from the triplets ii, jj, and vv such that L(ii(k),jj(k)) = vv(k).

if (nargout == 2)  && isWDefined
    DisR_Simulated = exp(-reshape(L*W(:), NTheta, []));
end

end



function [delta_d0, delta0] = getCORDrift(n_delta, ind_cr, WSz, NTheta, dxdy)
% set delta_d0, delta0: the center of rotaion (COR) drift: error, default no error, so set delta_d0, delta0 to all zero
driftfactor=4; % ??
cor=dxdy(1)*driftfactor/2*(WSz(1)/4-WSz(1)/3);  % XH: about -pixlSz*16.67
cr =[0 0; -cor -cor; -cor cor; cor cor; cor -cor;cor 0;0 cor;-cor 0;0 -cor];  % ??
cr = cr([1 4 6],:);  % XH: ?? why three cor shifts, probably not related to shift drift?

Delta_D = [0,1*cor];
delta_d0 = Delta_D(1);

if(n_delta == 2*NTheta)
    rng('default');
    pert = (-1+2*rand(NTheta,2))*dxdy(1)*2;
    pert = cumsum(pert,1)*2;
elseif(n_delta==4)
    pert1 = 0.5*dxdy(1)*10; pert2=2*dxdy(1)*10;
    pert = [pert1*ones(floor(NTheta/2),2);pert2*ones(NTheta-floor(NTheta/2),2)];
else % n_delta == 0
    pert = zeros(NTheta,2)*dxdy(1)*1;
end
if(ind_cr==2)
    delta0 = repmat(cr(ind_cr,:),NTheta,1)+pert*1;
    delta0(delta0==inf) = 0;
else
    delta0 = sparse(NTheta,2);
end
end

