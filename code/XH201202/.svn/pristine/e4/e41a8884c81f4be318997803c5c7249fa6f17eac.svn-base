function linkAxesXYZLimColorView(hf)
%set the same xyz limit/view/cameraPosition for all axes (such as created by subplot) in a figure
% Input
%   hf: handle of a figure, we will use current figure if hf is not specified
% Xiang Huang: xianghuang@gmail.com
% 2017.03.28: added to link the colormap 'Clim', 'CMode'
% 2018.12.21: added check property exist for all axes before link it. 
%   fixed the long time warning, this warning may crash the matlab, when close the image. 
%  ----------------
%  When call linkprop() in linkAxesXYZLimColorView()
%   Warning: Invalid property 
%   In localUpdateListeners (line 44)
%   In matlab.graphics.internal.LinkProp/set.Targets (line 72)
%   In matlab.graphics.internal.LinkProp (line 119)
%   In linkprop (line 33)
%   In linkAxesXYZLimColorView (line 27) 
%  ----------------
%  When close the figures, following warning and "MATLAB has encountered an internal problem and needs to close."
%   Warning: Invalid property 
%   In localUpdateListeners (line 44)
%   In matlab.graphics.internal.LinkProp/set.Targets (line 72)
%   In matlab.graphics.internal.LinkProp/removetargest (line 15)
%   In matlab.graphics.internal.LinkProp/processRemoveHandle (line 4)
%   In localUpdateListeners>@(varargin)hLink.processRemoveHandle(varargin{:}) (line 38)
%   In closereq (line 18)
%   In close
%   In close
% Todo: put this in Dropbox\Sourcecode_XHuang folder

if ~exist('hf', 'var') || isempty(hf)
    hf = gcf;
end

ax = findobj(hf, 'type', 'axes');
if numel(ax) < 2
    % for figure without axis or just one axis, no need to link, just return
    return;
end

%ax = hf.Children(arrayfun(@(x) isa(x, 'matlab.graphics.axis.Axes'), hf.Children)); % same as above, fin a figure could have children of (1) axes (2) UI objects, only select axes

% Find the axis limit that covers all axes. 
%Todo: verify this: for any axes, XLim, YLim, ZLim property always exist and has values?
%XLimAll: N*2 array, kth row has [xmin, xmax] for kth axes. It stacks ax(k).XLim which is a 1*2 array 
xLimAll = cat(1, ax(:).XLim); % Compare with  cat(2, ax(:).XLim), i.e., [ax(:).XLim], which is 1*2N array
yLimAll = cat(1, ax(:).YLim);
zLimAll = cat(1, ax(:).ZLim);

xLim = [min(xLimAll(:,1)), max(xLimAll(:,2))];
yLim = [min(yLimAll(:,1)), max(yLimAll(:,2))];
zLim = [min(zLimAll(:,1)), max(zLimAll(:,2))];

% Method I: linkprop, zoom in one axes will also automatically zoom in others (actually didn't work that way)h
% Link views so that rotation of view will remain in sync across both axes
%% a property is linked if it is the candidate and exist in all ax
propCandidates = {'XLim', 'YLim', 'ZLim', 'Color', 'View', 'CameraViewAngle', 'CLim', 'CMode'};
NCandidates = numel(propCandidates);
isExist = true(1, NCandidates);
for n = 1:numel(ax)
    % we always have ax >= 2, as ax < 2 will be returned above
    isExist = isExist & cellfun(@(x) isprop(ax(n), x), propCandidates);        
end
hLink = linkprop(ax, propCandidates(isExist));
%hLink = linkprop(ax, {'XLim', 'YLim', 'ZLim', 'Color', 'View', 'CameraViewAngle', 'CLim', 'CMode'}); % bad old code


%%
setappdata(hf, 'graphfics_linkaxes', hLink);
ax(1).XLim = xLim;
ax(1).YLim = yLim;
ax(1).ZLim = zLim;

% 2017.02.15_added set caxis the same as the left most sub figure
caxis0 = caxis(ax(end)); % the left most ax;
for n=1:numel(ax)
    caxis(ax(n), caxis0);
end

end
% Other Method II: simpler choice, but can't zoom in same time
% axis(ax, [XLim, YLim, ZLim]);