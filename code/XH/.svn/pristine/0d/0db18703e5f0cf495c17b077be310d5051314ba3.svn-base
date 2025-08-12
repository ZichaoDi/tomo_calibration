function [Theta, SourceKnot0, DetKnot0] = Define_Detector_Beam_Gaussian_XH(omega, dxdy, delta_d0, NTheta, NTau)
%Define detector and source beam
% omega = [-2     2    -2     2].*Tol; % XH: [xMin, xMax yMin yMax] in physical unit
% Output:
%   SourceKnot0:  source knot points, NTau*2 array, each row is one knot point's x,y coordinates
%   DetKnot0:     detectorlet knot points, NTau*2 array, each row is one knot point's x,y coordinates
%   Theta:        NTheta*1 scan angels in radians evenly distributed from 0 to pi-pi/NTheta. Note theta=0 and theta=pi are the same scan line 
%
% 2018.08.17:   XH refactored from script to function
% 2018.08.20:   XH refactored to eliminate duplcated inputs WSz and WSz(1)
% 2018.09.04:   XH replaced nTau+1 --> NTau (strange and easily confused with the +1), NTheta ---> NTheta
% 2018.11.30:   XH removed sythetic parameter, as it is not needed!


dTau = dxdy(1); % width of each discrete beam, or the gap between center of beamline. We assume the beamline gap is similar to object pixel size. Currently assume dx=dy, object pixel is square. 
halfDist = dTau*(NTau-1)/2; % half distance from source knot to destination knot
objDiag = sqrt(sum((omega([2,4])-omega([1,3])).^2)); % diagonal of the square object
assert(halfDist*2 >=  objDiag, 'Beamline should travel enough distance to cover the whole object, i.e., cover its diagonal.');

knot = dTau*(-(NTau-1)/2 : (NTau-1)/2)' + delta_d0; % NTau*1 column vector, plus the offset from delta_d0 amount drift of center of rotation (assume 0 currently)
  
SourceKnot0 = [repmat(+halfDist, NTau,1), knot];  % NTau*2 matrix, each row is the [x, y] coordinates
DetKnot0    = [repmat(-halfDist, NTau,1), knot];  % NTau*2 matrix, each row is the [x, y] coordinates


%%%=========== Assign Projection Angles;
%Theta = linspace(1, 360, NTheta) * pi/180;% XH: original code by Wendy, bad, as need 360 -> 180 to avoid duplicate, and start with 0 is more natural
Theta = (0:(NTheta-1)) * pi/NTheta; % same as linspace(0, 180-180/NTheta, NTheta);
