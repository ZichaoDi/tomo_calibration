function W = CreateElement_XH(WSz, NumElement, cross_ind)
%Create the object image, such as a 2D object image from central slice of 3D phantom
% Inputs:
% WSz:          size of W, 2D array [Ny, Nx].
% NumElement:   spectral bands, 1 at 2018.08.17
% cross_ind:    If not specified, we use the central slice of 3D phantom to generate 2D image.
% Output:
%
% 2018.08.17:   XH refactored from script to function
% 2018.08.20:   XH refactored to eliminate duplcated inputs WSz and WSz(1)

% A=double(imread('MRIhead-R.jpg'));
A = phantom3d(WSz(1));

if(~exist('cross_ind','var'))
    cross_ind = floor(WSz(1)/2);
end

A = squeeze(A(:,:,cross_ind));% sum(imread('phantom.png'),3);% XH: only use central 2D slice

% A=map1D(double(A),[0,1]);
if(NumElement==1)
    W=sum(A,3);% abs(peaks(WSz(1))); 
    % W(W<0.1)=0;
    % W(40:100,40:100)=W(40:100,40:100)+1e-2;
else
    % A=ones(WSz);
    A(A(:)<0)=0;
    tol=eps^(1/2);
    W=zeros(WSz(1),WSz(2),NumElement);
        val=unique(A(:));
    if(length(val)>1)
        % val=val(val~=0);
        val_i=[];
        for i_v=1:length(val)
            val_i(i_v)=length(find(A(:)==val(i_v)));
        end
        ExtraVal=[0.1 0.2 0.3];
        [~, i2]=sort(val_i,'descend');
        subind=[1 2 4];
        subind=subind(1:NumElement);
        if(length(val)<4)
            subind=[1 2];
            if(NumElement==1)
                subind=1;
            end
        end
        val=val(i2(subind));
        for i=1:length(subind)
            Ws=zeros(WSz(1),WSz(2));
            Ws(abs(A(:)-val(i))<tol)=val(i)+ExtraVal(i);
            W(:,:,i)=Ws;%+0.1;
        end
    end
end

end