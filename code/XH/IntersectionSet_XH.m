function [index,Lvec,linearInd] = IntersectionSet_XH(Source, Detector, xbox, ybox, theta,  ...
                                                   x, y, omega, WSz, dxdy, Tol, ...
                                                   plotTravel, BeforeEmit)
% Inputs:
%      plotTravel:  Optional, default false.  plot the intersection of beam with object
%      BeforeEmit:  Optional, default false.
% Functions called:
%      intersectLinePolygon()


finalfig = 99; fig2 = 99; fig5 = 99; % not sure the use of finalfig, fig2, fig5, though those were original set as global variable by Wendy, no other matlab files used those variables, maybe there were matlab files used those variables but not related to this tomography project. XH just set them to dummy values 99

if ~exist('plotTravel', 'var') || isempty(plotTravel)
    plotTravel = false; % Default false. NOT plot the intersection of beam with object
end
if ~exist('BeforeEmit', 'var') || isempty(BeforeEmit)
    BeforeEmit = false; % Default false
end



% [Ax, Ay] = polyxpoly([Source(1),Detector(1)],[Source(2),Detector(2)], xbox, ybox);
[intersects] = intersectLinePolygon([Source(1) Source(2) Detector(1)-Source(1) Detector(2)-Source(2)], [xbox',ybox']);
Ax=intersects(:,1);Ay=intersects(:,2);
if(isempty(Ax) || length(unique(Ax))==1 && length(unique(Ay))==1)
    % fprintf('no intersection \n')
    index=[];
    Lvec=[];
    linearInd=[];
else
    A=unique([Ax,Ay],'rows');Ax=A(:,1);Ay=A(:,2);
    if(theta==pi/2)
        Q=[repmat(Ax(1),size(y')),y'];
    elseif(theta==0 || theta==2*pi)
        Q=[x',repmat(Ay(1),size(x'))];
    elseif(theta==pi)
        Q=[x(end:-1:1)',repmat(Ay(1),size(x'))];

    elseif(theta==3*pi/2)
        Q=[repmat(Ax(1),size(y')),y(end:-1:1)'];
    else
        Q=[[x', (Ay(2)-Ay(1))/(Ax(2)-Ax(1)).*x'+(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1))];...
            [(y'-(Ay(1)*Ax(2)-Ay(2)*Ax(1))/(Ax(2)-Ax(1)))./((Ay(2)-Ay(1))/(Ax(2)-Ax(1))),y']];
    end

    indx=find(Q(:,1)-xbox(1)<-1e-6*Tol | Q(:,1)-xbox(3)>1e-6*Tol);
    indy=find(Q(:,2)-ybox(1)<-1e-6*Tol | Q(:,2)-ybox(2)>1e-6*Tol);
    Q=setdiff(Q,Q([indx;indy],:),'rows');
    Q=unique(Q,'rows');
    if(BeforeEmit)
        % dis=sqrt(sum(bsxfun(@minus,Q,Source).^2,2));
        % [~,InterOrder]=sort(dis);
        % Q=Q(InterOrder,:);
    else
        Q=unique([Q;A],'rows');
    end
    Lvec=sqrt(bsxfun(@minus,Q(2:end,1),Q(1:end-1,1)).^2 + bsxfun(@minus,Q(2:end,2),Q(1:end-1,2)).^2);
    %%%%%%%%%================================================================
    QC=(Q(1:end-1,:)+Q(2:end,:))/2;
    index=floor([(QC(:,1)-omega(1))/dxdy(1)+1, (QC(:,2)-omega(3))/dxdy(2)+1]);
    indInside=find(index(:,1)>0 & index(:,1)<=WSz(2) & index(:,2)<=WSz(1) & index(:,2)>0);
    index=index(indInside,:);
    [~,subInd]=unique(index,'rows');
    index=index(sort(subInd),:);
    Lvec=Lvec(indInside);
    Lvec=Lvec(sort(subInd));

    linearInd=sub2ind(WSz,index(:,2),index(:,1));
    %%%%%%%%%================================================================
    if plotTravel
        if(BeforeEmit)
            figure(finalfig)
            subplot(1,2,1);
            set(fig2,'visible','off');
            drawnow;
            if(~isempty(index))
                fig2 = plot((index(:,1)-1/2)*dxdy(1)-abs(omega(1)),(index(:,2)-1/2)*dxdy(2)-abs(omega(3)),'r*');%,Q(:,1),Q(:,2),'g-');
            end
        else
            figure(finalfig)
            subplot(1,2,1);
            set(fig5,'visible','off');
            drawnow;
            if(~isempty(index))
                fig5 = plot((index(:,1)-1/2)*dxdy(1)-abs(omega(1)),(index(:,2)-1/2)*dxdy(2)-abs(omega(3)),'bo',Q(:,1),Q(:,2),'g-');
            end
            pause;
        end
    end
    %%%%%%%%%================================================================

end
