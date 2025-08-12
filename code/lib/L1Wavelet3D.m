function [sparsity, sparsityNormalized] = L1Wavelet3D(x,optw)
% x: 3d matrix
% optw: choice of wavelet
% sparsity: scaler of the sparisty evaluation of x using wavelet defined by optw

wlev = optw.wlev;
wname = optw.wname;

wd = wavedec3(x,wlev,wname);
cw = vectorize_w(wd);

sparsity = sum(abs(cw(:)));
% normalize the sparsity with maximum absolute value and the number of elements
if nargout == 2
    t = abs(cw(:)); sparsityNormalized = sum(t/max(t))/numel(t);
end
end


function [cw] = vectorize_w(wd)

len_wd = length(wd.dec);
cw = []; 
for i = 1:len_wd
    cw = [cw; wd.dec{i}(:)];
end

end

