function [tbl_stats, NMIQR, p, X] = dotplottype_NG(Data, var, ylimit,ytitle)
    %Creates a dotplot of Data to compare up to 4 groups specified with T
    %   Axes can be specified to scale the annotations if the figure is
    %   part of a subplot
    %labels = {'Exc L2 LAMP5 LTK','Exc L2-3 LINC00507 FREM3','Exc L2-4 LINC00507 GLP2R','Exc L3-4 RORB CARM1P1','Exc L3-5 RORB COL22A1'} ;
    newlabel=[1 2 3 4 5 6];
    labels={'LTK', 'FREM3_s', 'FREM3_d', 'GLP2R','CARM1P1', 'COL22A1'};
    dat = Data(ismember(Data.NewLabel,newlabel),{'sample_id','NewLabel',var}) ;
    dat = sortrows(dat,'NewLabel') ;
    
    dat2 = table2array(dat(:,{var})) ;
    % select and prepare clusters
    [~,~,C] = unique(dat(:,2)) ;

    
    cl=unique(C(~isnan(C)));
    nclust=length(cl);
    for i=1:nclust
        T(C==cl(i))=i;
    end

    for i=1:nclust
        c(i)=length(T(T==i));
    end

    X=NaN(max(c),nclust);
    for i=1:nclust
        X(1:c(i),i)=dat2(T==i);
        %labels{i}=ClustNames{i};
        NMIQR(i,1) = sum(~isnan(X(:,i))) ;
        NMIQR(i,2) = nanmedian(X(:,i)) ;
        NMIQR(i,3) = quantile(X(:,i),0.25);
        NMIQR(i,4) = quantile(X(:,i),0.75);
    end
    
    if nclust == 2
        p = ranksum(X(:,1),X(:,2)) ;
    else
        [p,tbl,stats] = kruskalwallis(dat2,T,'off') ;
        %p = anova1(Data,T,'off');
    end
    
    realfig = figure('Position', [2000 10 1100 400]) ; %TDL 1800 width
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    
    subplot(1,2,1);
    if length(labels)==5
        dotdensity(X, 'medianLine', 'on', 'labels', labels, 'dotEdgeColor', [0.2 1 0.2; 0.5 0.8 1; 0.1 0.8 0.6; 0 0 1; 0.8 0.2 0.8], 'dotFaceColor', 'w','dotSize',5, 'lineColor','k' , 'lineWidth', 1.5, 'spread', 0.5);
    elseif length(labels)==6
        dotdensity(X, 'medianLine', 'on', 'labels', labels, 'dotEdgeColor', [0.2 1 0.2; 0.7 0.7 0.7; 0.5 0.8 1; 0.1 0.8 0.6; 0 0 1; 0.8 0.2 0.8], 'dotFaceColor', 'w','dotSize',7, 'lineColor','k' , 'lineWidth', 1.5, 'spread', 0.5);           
    elseif length(labels)==3
        dotdensity(X, 'medianLine', 'on', 'labels', labels, 'dotEdgeColor', [1 0.8 0.1; 0.8 0 0.15; 1 0.5 0.2], 'dotFaceColor', 'w','dotSize',5, 'lineColor','k' , 'lineWidth', 1.5, 'spread', 0.5);
    elseif length(labels)==8
        dotdensity(X, 'medianLine', 'on', 'labels', labels, 'dotEdgeColor', [0.2 1 0.2; 0.5 0.8 1; 0.1 0.8 0.6; 0 0 1; 0.8 0.2 0.8; 1 0.8 0.1; 0.8 0 0.15; 1 0.5 0.2], 'dotFaceColor', 'w','dotSize',5, 'lineColor','k' , 'lineWidth', 1.5, 'spread', 0.5);
    end  
    %dotdensity(X, 'medianLine', 'on', 'labels', labels, 'dotEdgeColor', [0.5 0.8 1; 1 0.4 0.4; 0.5 0.8 1; 1 0.4 0.4; 0.5 0.8 1; 1 0.4 0.4], 'dotFaceColor', 'w','dotSize',5, 'lineColor','k' , 'lineWidth', 1.5, 'spread', 0.5);
    %dotdensity(X, 'meanLine', 'on', 'labels', labels, 'dotEdgeColor', [0.7 0.7 0.7], 'dotFaceColor', 'w','dotSize',5, 'lineColor','k' , 'lineWidth', 1, 'spread', 0.5);
    set(gca, 'XLim', [0.5 max(T)+0.5]);
    if ~isempty(ylimit)
        ylim(ylimit)
    end
    box off
    if ~isempty(ytitle)
        ylabel(ytitle) ;
    else
        ylabel(var) ;
    end
    
    text (min(xlim)+diff(xlim)*0.05, min(ylim)+diff(ylim)*0.9, sprintf(['p: ' num2str(p,2)]), 'FontSize',14); 
    realax = subplot(1,2,2);

    tempfig = figure();     %make it current

    %%% !!!!
    % Statisisch toets
    
    comp = multcompare(stats) ;
    tbl_stats = array2table(comp,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
    tempax = tempfig.Children;
    set(tempax, 'Parent', realfig, 'Position', realax.Position);  %move it between figures
    delete(realax);
    close
    
  
    %savename = ['Dot_' var] ;
    savename = ['Dot_Split_' var] ;
%    print(fullfile('C:\Users\DBHeyer\Documents\PhD\RNA analysis\Patchseq\Manuscript',[savename '.eps']),  '-depsc','-r300','-painters');
    %print(fullfile('C:\Users\DBHeyer\Documents\PhD\RNA analysis\Patchseq\Figures\genesets\Dotplots',[savename '.jpg']),  '-djpeg', '-r300');   
    %savefig(fullfile('C:\Users\DBHeyer\Documents\PhD\RNA analysis\Patchseq\Figures\genesets\Dotplots\Interactive',[savename '.fig']));   
end


function [varargout] = dotdensity(varargin)
%dotdensity   Dot density plot
%   dotdensity(X,Y) plots a dot density plot of the columns of Y versus X.
%   Each column of y is plotted centered on the corresponding value of X.
%   
%
%   dotdensity(Y) plots the columns of Y versus their index.
%
%   dotdensity(AX,....) plots into the axis with the handles AX.
%
%   dotdensity(..., 'PARAM1', val1, 'PARAM2', val2....) specifies
%   optional parameter/value pairs.
%       'medianLine'  'on' draws a horizontal line at the median of each 
%                      data set
%
%       'meanLine'     'on' draws a horizontal line at the mean  of each 
%                      data set.
%
%       'labels'       Cell array of strings, or numeric vector, of labels  
%                      for each data set.
%                   
%       'delta'        A scalar (default is 0.02) that sets the spacing
%                      between x-coordinates of one column of Y.
%
%       'spread'       A scalar (default is 0.49 * distance between
%                      the first two X values) indicating the maximum 
%                      spread of one column of Y from the corresponding X
%                      value
%
%       'nBins'        An integer (default is 10). A histogram of NBINS is 
%                      used to divided yCol in groups. Points in the same 
%                      bin are assumed to be likely to overlap in the
%                      y-direction and are therefore distributed in the 
%                      x-direction with a distance 'delta' and a maximum
%                      spread 'spread.'
%
%       'dotEdgeColor' Specifies edge color of markers used for dots. This
%                      value can be a string ('g') or three column numeric 
%                      array ([0.1,0.2,0.3]) of RGB values and can 
%                      specify single colors or multiple.  If multiple
%                      colors are specified each column of Y is drawn in
%                      a different color.
%
%       'dotFaceColor' Specifies the face color of markers used for dots. 
%                      Uses the same format as 'dotEdgeColor'  
% 
%       'dotMarker'    A string that specifies the marker used for dots.
%                      The default is 'o'. 
%
%       'dotSize'     The size of the marker used for dots.
% 
%       'lineColor'    The color of the median or mean line. The format
%                      either a string ('g') or a three element vector 
%                      ([1,0,0]) is accepted. The default is 'b'.
%
%       'lineWidth'    The width of the median or mean line. 
%
%   HDots = dotdensity(..) returns a vector of handles to dot density  
%   plots.
%
%   [HDots, HLines] = dotdensity('Median','on') or 
%   [HDots, HLines] = dotdensity('Mean','on') returns a vector of handles      
%   to the mean or median line.
%
%   Copyright (c) Molly J. Rossow <molly@mollyrossow.com> 8/1/2013
%
%   Revision history
%   12/4/2013 Fixed bug in findAx so the first argument to dotdensity 
%             can be a one element vector. 
% 
%   Example: 
%   x = [10, 20, 30, 40];
%   m = [4, 2, 8, 5];
%   sd = [0.1, 0.5, 0.7, 2];
% 
%   y = zeros(100, 4);
%   for c = 1:length(x)
%       y(:,c) = m(c) + sd(c) * randn(100, 1);
%   end
% 
%   figure
%   dotdensity(x, y, 'medianLine', 'on')
% 


%%
% Parse arguments.
[x, y, ax, optionalArgs] = parseArgs(varargin);

% Set spread relative to x-coordinates.
if optionalArgs.spread == 0
    if length(x) == 1
        spread = 0.4;
    else
        spread = 0.4 * abs(x(1) - x(2));
    end
else
    spread = optionalArgs.spread;
end

% Set delta relative to x-coordinates.
if optionalArgs.delta == 0
    if length(x) == 1
        delta = 0.02;
    else
        delta = spread * 0.1;
    end
else
    delta = optionalArgs.delta;
end

%%
% Calculate new x coordinates and rearrange y coordinates.
newX = zeros(size(y));
newY = zeros(size(y));
for ii = 1:length(x)
    yCol = y(:,ii);
    [newX(:,ii), newY(:,ii)] =generateXCoords(x(ii),yCol,...
        delta, spread, optionalArgs.nBins);
end

%%
% Draw plot

% Draw dot density plot.
dotHandles = plot(ax, newX, newY, 'o','LineWidth',1.5);

% Set edge color of dots.
if ~isempty(optionalArgs.dotEdgeColor)
    setEdgeColor(optionalArgs.dotEdgeColor, dotHandles)
end

% Set face color of dots.
if ~isempty(optionalArgs.dotFaceColor)
    setFaceColor(optionalArgs.dotFaceColor, dotHandles)
end

% Set marker for dots.
if ~ isempty(optionalArgs.dotMarker)
    setMarker(optionalArgs.dotMarker, dotHandles)
end

% Set dot size.
if optionalArgs.dotSize ~= 0
    setDotSize(optionalArgs.dotSize, dotHandles)
end


% Draw median line.
lineHandles = [];
if strcmp('on', optionalArgs.medianLine)
    h = drawMedianLines(x,y,spread);
    lineHandles = [lineHandles, h];
end

% draw mean line
if strcmp('on', optionalArgs.meanLine)
    h = drawMeanLines(x,y,spread);
    % Concatenate h with lineHandles so all handles lines are returned
    % together.
    lineHandles = [lineHandles; h];  
end

% Set line color for median and/or mean line.
if (~isempty(optionalArgs.lineColor) && ~isempty(lineHandles))
    setLineColor(optionalArgs.lineColor, lineHandles)
end

% Set line width for median and/or mean line.
if (~isempty(optionalArgs.lineWidth) && ~isempty(lineHandles))
    setLineWith(optionalArgs.lineWidth, lineHandles)
end

% Add lables.
if ~isempty(optionalArgs.labels)
    set(ax,'XTick',x)
    set(ax,'XTickLabel',optionalArgs.labels)
    %set(ax,'XTickLabel', [])
    xtickangle(30)
end

% Define varargout.
if nargout == 1
    varargout{1} = dotHandles;
elseif nargout == 2
    varargout{1} = dotHandles;
    varargout{2} = lineHandles;
end
end

%% 
% Set line width for mean and/or median lines.
function setLineWith(lineWidth, handles)
set(handles,'LineWidth',lineWidth)
end

%% 
%  Set line color for mean and/or median lines.
function setLineColor(lineColor, handles)
   index = 1;
if ischar(lineColor)
    for h = handles'
        set(h,'Color',lineColor(index))
        % Increment index to cycle through colors.
        if index < length(lineColor)
            index = index + 1;
        else
            index = 1;
        end
    end
elseif isnumeric(lineColor)
    s = size(lineColor);
    if s(2) ~= 3 % if LineColor is not a string or 3 column array
        error('LineColor must be a three column array or a string')
    else
        for h = handles'
            set(h,'Color',lineColor(index,:))
            if index < s(1)
                index = index + 1;
            else
                index = 1;
            end
        end
    end
end
end 

%% 
% Set the dot size.
function setDotSize(size, handles)
set(handles, 'MarkerSize', size)
end

%% 
% Set the dot marker,
function setMarker(markers, handles)
index = 1;
% Set each handle seperately in case there are multiple markers.
for h = handles'
    set(h, 'Marker', markers(index))
    % Increment index to cycle through markers.
    if index < length(markers)
        index = index + 1;
    else
        index = 1;
    end
end
end

%% 
% Set dot edge colors.
function setEdgeColor(edgeColor, handles)
index = 1;
if ischar(edgeColor)
    for h = handles'
        set(h,'MarkerEdgeColor',edgeColor(index))
        % Increment index to cycle through colors
        if index < length(edgeColor)
            index = index + 1;
        else
            index = 1;
        end
    end
elseif isnumeric(edgeColor)
    s = size(edgeColor);
    if s(2) ~= 3 % if edgeColor is not a string or 3 column array
        error('DotEdgeColor must be a three column array')
    else
        for h = handles'
            set(h,'MarkerEdgeColor',edgeColor(index,:))
            if index < s(1)
                index = index + 1;
            else
                index = 1;
            end
        end
    end
end
end

%% 
% Set dot face color.
function setFaceColor(faceColor, handles)
index = 1;
if ischar(faceColor)
    for h = handles'
        set(h,'MarkerFaceColor',faceColor(index))
        % Increment index to cycle through colors.
        if index < length(faceColor)
            index = index + 1;
        else
            index = 1;
        end
    end
elseif isnumeric(faceColor)
    s = size(faceColor);
    if s(2) ~= 3 % if edgeColor is not a string or 3 column array
        error('dotFaceColor must be a three column array')
    else
        for h = handles'
            set(h,'MarkerFaceColor',faceColor(index,:))
            if index < s(1)
                index = index + 1;
            else
                index = 1;
            end
        end
    end
end
end


%%
% Draw median lines for each column of Y.
function h = drawMedianLines(x,y,spread)
h = zeros(length(x),4);
% Reshape to be a column vector.
% h = reshape(h, [length(h), 1]);
for ii = 1:length(x)
    yCol = y(:,ii);
    m = median(yCol,'omitnan'); 
    % IQR
    lqrt = quantile(yCol,0.25);
    hqrt = quantile(yCol,0.75);  
    h(ii,1) = line([x(ii) - spread*0.6, x(ii) + spread*0.6],[m,m]);
    h(ii,2) = line([x(ii),x(ii)],[lqrt,hqrt]);
    h(ii,3) = line([x(ii) - spread*0.3, x(ii) + spread*0.3],[hqrt,hqrt]);
    h(ii,4) = line([x(ii) - spread*0.3, x(ii) + spread*0.3],[lqrt,lqrt]);
    % SD
%     sd = std(yCol,'omitnan');
%     h(ii,1) = line([x(ii) - spread*0.6, x(ii) + spread*0.6],[m,m]);
%     h(ii,2) = line([x(ii),x(ii)],[m-sd,m+sd]);
%     h(ii,3) = line([x(ii) - spread*0.3, x(ii) + spread*0.3],[m+sd,m+sd]);
%     h(ii,4) = line([x(ii) - spread*0.3, x(ii) + spread*0.3],[m-sd,m-sd]);
end
end

%% 
% Draw mean lines for each column of Y.
function h = drawMeanLines(x,y,spread)
h = zeros(length(x),4);
% h = reshape(h, [length(h), 1]); removed after SD lines addition (R.W.)
for ii = 1:length(x)
    yCol = y(:,ii);
    m = mean(yCol,'omitnan'); 
    sd = std(yCol,'omitnan');
    h(ii,1) = line([x(ii) - spread*0.6, x(ii) + spread*0.6],[m,m]);
    h(ii,2) = line([x(ii),x(ii)],[m-sd,m+sd]);
    h(ii,3) = line([x(ii) - spread*0.3, x(ii) + spread*0.3],[m+sd,m+sd]);
    h(ii,4) = line([x(ii) - spread*0.3, x(ii) + spread*0.3],[m-sd,m-sd]);
end
end

%% 
% Parse the input functions including optional arguments. 
function [x, y, ax, optional] = parseArgs(inArgs)
% Select default values for optional arguments. 
defaultMedianLine = 'off';
defaultMeanLine = 'off';
defaultLabels = {};
defaultDelta  = 0;
defaultSpread = 0;
defaultNBins = 10;
defaultEdgeColor = '';
defaultFaceColor = ''; 
defaultMarker = 'o';
defaultDotSize = 0; 
defaultLineColor = '';
defaultLineWidth = ''; 

% If the following functions find the argument they are looking for they
% remove it from inArgs
[ax, inArgs] = findAx(inArgs); % Returns the argument ax or the current 
                               % axis.
[x, y, inArgs] = findXandY(inArgs); % x is the input argument X if it 
                                    % exists and a vector of indices 
                                    % otherwise.

%%
% Parse optional parameter value pairs. 
p = inputParser;
p.CaseSensitive = false;
addParamValue(p,'medianLine',defaultMedianLine,@ischar);
addParamValue(p,'meanLine',defaultMeanLine,@ischar);
addParamValue(p,'delta',defaultDelta, @isscalar);
addParamValue(p,'spread',defaultSpread, @isscalar);
addParamValue(p,'nBins',defaultNBins, @isscalar);
addParamValue(p,'labels',defaultLabels); 
addParamValue(p, 'dotEdgeColor', defaultEdgeColor);
addParamValue(p, 'dotFaceColor', defaultFaceColor);
addParamValue(p, 'dotMarker',defaultMarker);
addParamValue(p, 'dotSize', defaultDotSize);
addParamValue(p, 'lineColor', defaultLineColor);
addParamValue(p, 'lineWidth', defaultLineWidth);

parse(p,inArgs{:});
optional = p.Results; 
end

%%
% Extract optional argument ax. Ax must be the first argument in INARGS.
function [ax, inArgs] = findAx(inArgs)
ax = gca;
if isnumeric(inArgs{1})
    if numel(inArgs{1}) == 1
        if ishandle(inArgs{1})
            if strcmp('axes',get(inArgs{1},'type'))
                ax = inArgs{1};
                % Remove ax from inArgs so it can be ignored when parsing 
                % other arguments.
                inArgs = inArgs(2:end); 
            end
        end
    end
end
end

%%
% Check first and second input arguments (after ax has been removed, if it
% was there) to see if they are appropriate values for X and Y. 
function [x, y, inArgs] = findXandY(inArgs)

% Extract mandatory argument Y and optional argument X
% Y might be vector or matrix and can be the first or second input
if length(inArgs) >= 2
in1 = inArgs{1};
in2 = inArgs{2};
if isnumeric(in1) && isnumeric(in2) % in1 and in2 are two numerical arrays,
                                    % X and Y.
    if ~isvector(in1) % If x is not a vector.
        error(['X is of size ', num2str(size(in1)), 'X must be a vector.']);
    else % If x is a vector
        lx = length(in1);
        sy = size(in2);
        % If the # of columns of Y and elemnts of x don't match.
        if ~(lx == sy(2)) 
            error(['The number of columns in Y must be the same as the '...
            'number of values in X']);
        else % If X and Y are valid arrays for plotting. 
            x = in1;
            y = in2;
            % Remove X and Y from inArgs.
            if length(inArgs) > 2
                inArgs = inArgs(3:end);
            elseif length(inArgs) == 2
                inArgs = {};    
            end
        end
    end
else % X and Y are not both present or valid arguments. 
    if isnumeric(in1)
        y = in1;
        % Define x as indices of Y
        s = size(y);
        x = 1:1:s(2);
        % remove Y from inArgs
        if length(inArgs) > 1
            inArgs = inArgs(2:end);
        elseif length(inArgs) == 1
            inArgs = {};
        end
    else
         error('Must have numeric input for Y')
    end
end
elseif length(inArgs) == 1 % Only Y is given as an input argument. 
 in1 = inArgs{1};
  % Define x as indices of Y
  if isnumeric(in1)
      y = in1; 
        s = size(y);
        x = 1:1:s(2);
      inArgs = {};
  else
      error('Must have numeric input for Y')
  end
else
    error('Must have numeric input for Y')
end
end

%% 
% Creates a a set of x-coordinates for plotting one column of a dot density
% plot. Returns NEWX a vector of x-coordinates centered on X with spacing
% DELTA for plotting with NEWY. NEWY is a reordered version YCOL. Values in
% newX are never more than MAXSPACE from X to prevent overlap with other
% plots. A histogram of NBINS is used to divide yCol in groups. Points in
% the same bin are assumed to be likely to overlap in the y-direction and
% are therefor distributed in the x-direction.
function [newX, newY] = generateXCoords(x,yCol,delta,maxSpace,nBins)

% Remove nan values from yCol. Keep track of the original size of yCol.
sYCol = size(yCol);
yCol = yCol(~isnan(yCol));

% Initialize newX and newY.
newX = zeros(size(yCol));
newY = zeros(size(yCol));

% Calculate the maximum number of elements in row of spread x-values.
maxElementsInRow = floor(maxSpace / delta);

% Calculate histogram of yCol and sort yCol so y-coordinates can be
% associated with the appropriate x-coordiantes.
h = hist(yCol,nBins);
yCol = sort(yCol);

index = 1; % Index of elements added to returned vectors.

% Loop over histogram bins organizing y-values and creating x-values for
% each.
for numElements = h
    %Get the y-coordinates associated with this bin.
    tempY = yCol(index:index + numElements -1); 
    % Put them in a random order.
    irand = randperm(numElements); 
    % Add them to the output vector.
    newY(index:index + numElements - 1) = tempY(irand);
    % Create corresponding x-coordinates.
    if numElements == 1
        newX(index) = x;
    elseif numElements > 1
        tempX = createXVec(x, numElements, delta, maxElementsInRow);
         % Add x-coordinates to output vector
        newX(index: index + numElements - 1) = tempX;
    end 
    index = index + numElements;
end

% Pad yCol and x with nan to restore ycol to it's original size.
newSize = size(newY);
newY = [newY; nan*ones(sYCol(1) - newSize(1),1)];
newX = [newX; nan*ones(sYCol(1) - newSize(1),1)];
end

%% 
% Returns a vector of x-coordinates for all y-coordinates in one bin of a
% histogram. X-coordinates are centered on X, have a spacing of DELTA or 0
% and are no more than MAXELEMENTS*DELTA from the center value. 
function outVec = createXVec(x, num, delta, maxElements)  
nRows = ceil(num / maxElements); % Data plotted with these x-coordinates
                                 % willappear approximately in rows since
                                 % the x-coordinates are repeated (if
                                 % maxElements is large enough)and the y
                                 % values are close together.
outVec = [];
lastRowElements =  rem(num,maxElements);
if lastRowElements == 0
    lastRowElements = maxElements;
end
for ii = 1:nRows
    if ii == nRows % last row, fill with all remaining elements. 
        outVec = [outVec, createSingleXRange(x, lastRowElements, delta)];
    else % previous rows
        outVec = [outVec, createSingleXRange(x, maxElements, delta)];
    end
end
end

%% 
% Returns a vector of NUM elements centered on X with spacing Delta. The
% vector returned by this function corresponds to a single, approximate row
% of spread x-values in the final plot.
function outVec = createSingleXRange(x, num, delta)
outVec = x-delta*floor((num/2)):delta:x+delta*floor((num/2));
if mod(num,2) == 0 % if num is even
    outVec = removeMiddleElement(outVec);
end
end

%% 
% Removes the middle element of an odd length vector.
function outVec = removeMiddleElement(inVec)
midPoints = ceil(length(inVec) /2);
outVec = [inVec(1:midPoints - 1), inVec(midPoints + 1:end)];
end

