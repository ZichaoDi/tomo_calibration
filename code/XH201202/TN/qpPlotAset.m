function [] = qpPlotAset(Aset,iter,n,changeSize,nitOld)

%-----------------------------------------------------------------------
% [] = qpPlotAset(Aset,iter,n)
% 
% INPUT ARGUMENT:
% Aset = +1 lower, -1 upper bound
% iter = iteration number
% n    = number pf variables
%
% AUTHORS
%    Michael P. Friedlander     and     Sven Leyffer
%    michael@mcs.anl.gov                leyffer@mcs.anl.gov
%                  Argonne National Laboratory
%
% DEVELOPMENT
% 01 Jun 2004: First version derived from PlotAset.m.
%
% $Revision: 1.1.1.1 $ $Date: 2005/10/27 16:47:48 $
%-----------------------------------------------------------------------

  % ... initialize/reset active c/s picture
%   if (iter == 1) 
%      clf reset; 
     title('Cauchy Active Sets (red=low, grn=up, )','FontSize',10);
     ylabel('variable index','FontSize',8); 
     xlabel('iteration','FontSize',8);
     axis([0 10+iter 0 n]); hold on;
%   elseif (mod(iter,10) == 0)
%      hold off; axis([0 iter+10 0 n]); hold on;
%   end;  

  % ... plot active constraints
  for i=1:length(Aset)
     ix = [ iter-0.5 iter+0.5 iter+0.5 iter-0.5 ];
     iy = [ i-0.5    i-0.5    i+0.5    i+0.5 ];
     if (Aset(i) == -1)
	 fill(ix,iy,'r','LineStyle','none'); 
     elseif (Aset(i) == 1)
	 fill(ix,iy,'g','LineStyle','none');     
     end;
     hold on;
  end;
 text((iter+nitOld)/2,length(Aset)/2,num2str(changeSize'),'FontSize',12);
  
      %%%#####################################visulize active set
    %     figure(44);
    %     subplot(1,2,1);
    %     qpPlotAset(ipivot,nit,length(x));
    %     [xBinding]=crash (x_new, low, up);
    %     subplot(1,2,2);
    %     qpPlotAset(xBinding,nit,length(x));
    %     drawnow;
    %      for ib=1:length(xBinding)
    %          if(xBinding(ib)==-1)
    %          MarkerStr{ib}='v';
    %          MarkerC{ib}='r';
    %          elseif(xBinding(ib)==1)
    %              MarkerStr{ib}='^';
    %              MarkerC{ib}='y';
    %          else
    %              MarkerStr{ib}='.';
    %              MarkerC{ib}='b';
    %          end
    % %          if(nit==0)
    % %          Track=figure('name','Track');
    % %          end
    %          figure(14);
    %          plot(nit,ib,'Marker',MarkerStr{ib},'MarkerFaceColor',MarkerC{ib});
    %          hold on;
    %          drawnow;
    %      end
    %%%##############################################################
