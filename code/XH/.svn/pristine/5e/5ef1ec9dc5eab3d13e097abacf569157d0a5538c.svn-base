function multAxes(funHandle, varargin)
%multAxes apply funHandle (@title/@imshow/@imagesc/@colorbar/axis etc.) to mutliple axes of current figure.
%It will create a new figure (if no figure exist) or new axes (if current figure has no axes)
%
% Input:
%   funHandle: such as @title/@imshow/@imagesc/@colorba
%   varargin:  input arguments for funHandle, each axes could have one
% Output:
%   NONE, only show figures
%
% Usage:
%   multAxes(@colorbar);  % applies color bar to all the axes in current figure
%   multAxes(@axis, 'image');  % set axes('image') for all aces in current figure
%   multAxes(@imshow, I1, I2);   or  multAxes(@imshow, {I1, I2});  % shows two images in current figure, e.g., left/right images
%   multAxes(@title, 'Left', 'Right');   or   multAxes(@title, {'Left', 'Right'}); % add title to multiple axes in current figure
%   multAxes(@imshow, {I1, I2}, {[Low1, High1], [Low2, High2]}) ; % show two images with different display range
%   multAxes(@imshow, {I1, I2}, [Low, High]);  % show two images with the same display range
%   multAxes(@imshow, {I1, I2}, [Low, High]); multAxes(@colorbar); linkAxesXYZLimColorView(); 

% Usage Explain: 
% (0) When funHandle takes no argument
%   multAxes(@colorbar);
%
% (1) When funHandle only takes one input argument
% (1.1) varargin{n} is the input argument for nth axes, when varargin is 1*N cell array
%   multAxes(@imshow, I1, I2);
%   multAxes(@title, 'Left', 'Right');
% (1.2) varargin{1}{n} is the input argument for nth axes, when varargin is 1*1 cell
%   multAxes(@imshow, {I1, I2});
%   multAxes(@title, {'Left', 'Right'});
%
% (2) When funHandle takes K (K>1) input arguments, varargin is 1*K cell and varargin{k} is kth input argument
% (2.1) varargin{k}{n} is the kth input argument for nth axes, when varargin{k} is 1*N cell array (e.g. {I1, I2}), 
%   multAxes(@imshow, {I1, I2}, {[Low1, High1], [Low2, High2]}); % show two imagew with different display range
% (2.2) varargin{k} is the kth input argument and the same for all the axes
%   multAxes(@imshow, {I1, I2}, [Low, High]); % show two images with the same display range
% NOTE: we don't recommend use multAxes(@imshow, {I1, I2}, {[Low, High]});
%
%
% References: arrayfun, cellfun, 
%           multAxes will replace titleMult/imshowMult/plotMult functions
%
% For example below works for title as title(AX, ...) correct but doesn't work for imshow etc as imshow(ax, ...) incorrect
% arrayfun(@title, axArray(:), varargin(:)); % 
%  
% Author: Xiang Huang: xianghuang@gmail.com
% 2017.02.21
%
if nargin == 0; testMe(); return; end % Unit Test


% (1) Setup figure: if there is no figure, create an figure, otherwise we use the exist current figure
if isempty(allchild(0))
    figure();
end

% (2) Findout N and K from varargin: total number of axes and function arguments.
% e.g. varargin =  {{I1,I2,I3}, [Low1, High1]} tells N = 3, K=2
N = 1;  % Start with unkown number of axes unless varargin tells the number of axes (>1)
if isempty(varargin)    
    % e.g. multAxes(@colorbar);
    K = 0;      
elseif numel(varargin) == 1  
    % e.g. multAxes(@imshow, {I1, I2}) or multAxes(@imshow, I);
    K = 1;  
    if iscell(varargin{1})
        N = numel(varargin{1});        
    end
elseif ~any(cellfun(@iscell, varargin))   
    % e.g. multAxes(@imshow, I1, I2);
    K = 1;
    N = numel(varargin);
    varargin = {varargin}; % convert multAxes(@imshow, I1, I2) to standard form of multAxes(@imshow, {I1, I2});
else
    % e.g. multAxes(@imshow, {I1, I2}, [Low, High]); 
    K = numel(varargin); % K > 1 always
    for i = 1:numel(varargin)
        if iscell(varargin{i}) && numel(varargin{i}) ~= 1  % only process varargin{i}, if it is {I1, I2, ...} not {I} or I
            if N == 1
                N = numel(varargin{i});
            else
                assert(N==numel(varargin{i}), sprintf('Input arguments indicates both %d and %d axes, inconsist', N, numel(varargin{i})));
            end            
        end
    end    
end


% (3) Set up axes
axArray = findobj(gcf, 'type', 'axes', '-depth', 1); % only check one depth
if (~isempty(axArray)) && (N==1 || N==numel(axArray))
    % (3.1) Draw on current axes, if already exist axes & current number of axes matches argument size
    % reverse the order so that the first one is the first axes        
    axArray = axArray(end:-1:1);  % reverse the order so that the first one is the first axes, but this only applies to 1D grid        
else
    % (3.2) Create N axes according to the number of input, if NO axes or current number of axes doesn't match argment size
    if (~isempty(axArray))  % && (N>1 && N~=numel(axArray))
        % Warning if overwrite current axes.
        warning('Over-write current axes: as input arguments indicates %d axes which doesnot match current figure with %d axes', N, numel(axArray));
    end

    % (3.2.1) determine the grid for sub-figures
    R = []; C = [];
    % (3.2.1.1) R*C grid if any input argument varargin{k} is 2D cell array of size R*C
    % If multiple input argument varargin{k} are 2D cell array, use the size of last one
    for k = 1:K
        if iscell(varargin{k}) && ~isvector(varargin{k})
            assert(ndims(varargin{k}) <= 2, sprintf('Each input argument for multiple axes need to be packed either in 1D or 2D array, currently %dD!', ndims(varargin{k})))
            [R, C] = size(varargin{k});
            varargin{k} = varargin{k}'; % Transpose as subplot(R,C,n) is row-major order, but array f/voxSz/th are column-major order.
        end
    end
    % (3.2.1.2) If all input arguments are 1D vector, decide grid size from screen size
    if isempty(R) || isempty(C)
        allMonitorPositions = get(0,'MonitorPositions'); % could be one or mutliple monitors, a row for each monitor with [x0, y0, width, height]
        monitorAspectRatio = allMonitorPositions(1,3) / allMonitorPositions(1,4); % main/1st monitor width vs height
        if monitorAspectRatio < 1
            % vertical screen, user may use a widescreen by rotate it by 90 degree
            C = round(sqrt(N));
            R = ceil(N/C);
        else
            % widscreen, typical setup
            R = round(sqrt(N));
            C = ceil(N/R);
        end
    end
    % (2.2) create axes
    axArray = gobjects(N, 1);
    for n = 1:N
        axArray(n) = subplot(R,C,n);
    end
end 


% (4) apply funHandle for each axes
for n = 1:numel(axArray)
    axes(axArray(n));  % works for title/image/imagesc etc, but maybe slower if use title(axArray(n), varargin{n})
    if K == 0
        funHandle(); % colorbar() etc where no input argument
    else
        varList = cell(1, K);
        for k = 1:K
            if iscell(varargin{k}) 
                if numel(varargin{k})>1
                    % e.g. varargin{k} = {I1, I2} from multAxes(@imshow, {I1, I2}, [Low, High])
                    varList{k} = varargin{k}{n};
                else
                    % e.g. varargin{k} = {[Low, High]}; in case used (not recommended) multAxes(@imshow, {I1, I2}, {[Low, High]})
                    varList{k} = varargin{k}{1};
                end
            else
                % e.g. varargin{k} = [Low, High];  from multAxes(@imshow, {I1, I2}, [Low, High])
                varList{k} = varargin{k};                
            end
            
        end
        funHandle(varList{:});
    end
end

% linkAxesXYZLimColorView();

end


%% Unit Test
function testMe()
%%
I1 = imread('cameraman.tif'); I1 = double(I1); I1 = I1/max(I1(:)); 
I2 = fliplr(I1);
%%
set(figure(1), 'Name', 'Test with imshow and title using two axes');  clf; 
multAxes(@imshow, {I1, I2});
multAxes(@title, 'Left', 'Right');
%%
set(figure(2), 'Name', 'Test with imshow and title using six axes with user specified 3x2 grid'); clf;
multAxes(@imshow, {I1, I2; I1, I2; I1, I2});  % 3*2 subimages
multAxes(@title, '1', '2', '3', '4', '5', '6');

%%
set(figure(3), 'Name', 'Test with imagesc which takes two inputs, and colorbar which takes zero inputs'); clf;
multAxes(@imagesc, {I1, I2}, {[0, 1], [0, 0.5]}); multAxes(@colorbar); multAxes(@axis, 'image');
%%
set(figure(4), 'Name', 'Test with imagesc which takes two inputs, and one argument is the same for all axes'); clf;
multAxes(@imagesc, {I1, I2, I1}, [0, 0.5]); multAxes(@colorbar); multAxes(@axis, 'image');
%%
set(figure(5), 'Name', 'Test with imagesc which takes two inputs, and one argument warpped in cell(1) is the same for all axes'); clf;
multAxes(@imagesc, {I1, I2, I1}, {[0, 0.5]}); multAxes(@colorbar); multAxes(@axis, 'image');
%%
set(figure(6), 'Name', 'Test with imagesc for two images with different display range, then use linkAexsXYZLimColorView() to set display range same as left image'); clf;
multAxes(@imagesc, {I1, I2}, {[0, 1], [0, 0.5]}); multAxes(@colorbar); linkAxesXYZLimColorView();
end