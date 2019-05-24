function [h, hb, hscat] = scatterBox(x,grps,varargin)
%[h, hb, hscat] = scatterBox(x,grps,varargin)
%
%   I usually like to plot individual data points on top of boxplots so
%   this function does that by separating your x vector into grps and then
%   jittering the data points to overlay. 
%
%   INPUTS
%       x: vector of data points to plot.
%
%       grps: grouping vector the same size as x. Each value indicates
%       which group each value of x belongs to. 
%
%       NAME,VALUE arguments:
%           xLabels: cell array of strings specifying the group names, must
%           be in the same order as grps (i.e., the min of grps must be the
%           first string in xLabels and the max of grps must be the last. 
%           
%           yLabel: string specifying y axis label. 
%
%           boxColor: string or RGB value specifying what color you want
%           the boxplot to be. default = 'k' (black)
%
%           circleSize: scalar or vector specifying scatter plot circle
%           sizes. default = 20;
%
%           circleColors: string or RGB values specifying what color you
%           want all or individual circles in the scatter plot to be.
%           default = [0.7 0.7 0.7];
%       
%           transparency: scalar specifying how transparent you want the
%           circles to be. 1 is opaque, 0 is completely transparent.
%           Default = 0.3.
%
%           sf: spread factor, scalar specifying how wide you want the
%           individual data points to be jittered. Default = 0.05;
%
%           position: vector specifying the position of the figure. 
%
%           plotBox: whether or not to make boxplot.
%
%           h: axes to plot into (default = make new figure)
%
%           paired_id: if you are using paired data points, will connect
%           dots between each point with the same id # in each group.
%
%   OUTPUTS
%       h, hb, hscat: handles to axes, boxplot, and scatterplot
%       respectively
%

%% Set up. 
    p = inputParser; 
    p.addRequired('x',@(x) isnumeric(x));
    p.addRequired('grps',@(x) isnumeric(x));
    p.addParameter('xLabels',arrayfun(@(a) ['Group ' num2str(a)], unique(grps),...
        'UniformOutput',false),@(x) iscell(x));
    p.addParameter('yLabel','Metric',@(x) ischar(x)); 
    p.addParameter('boxColor','k',@(x) ischar(x) || isnumeric(x));
    p.addParameter('circleSize',20,@(x) isnumeric(x)); 
    p.addParameter('circleColors',[0.3 0.3 0.3],@(x) ischar(x) || isnumeric(x));
    p.addParameter('transparency', 0.7, @(x) isscalar(x)); 
    p.addParameter('sf',.05,@(x) isscalar(x));
    p.addParameter('position',[520 350 300 450]); 
    p.addParameter('plotBox',true,@(x) islogical(x));
    p.addParameter('paired_ind', [], @(a) isempty(a) || isnumeric(a));
    p.addParameter('h', [], @(a) ishandle(a) || isempty(a)); 
    
    p.parse(x,grps,varargin{:});
    xLabels = p.Results.xLabels; 
    yLabel = p.Results.yLabel; 
    boxColor = p.Results.boxColor;
    circleSize = p.Results.circleSize;
    transparency = p.Results.transparency;
    circleColors = p.Results.circleColors;
    sf = p.Results.sf; 
    position = p.Results.position;
    plotBox = p.Results.plotBox;
    h = p.Results.h;
    paired_ind = p.Results.paired_ind; 
    
    %Turn into column.
    if size(x,1) < size(x,2) 
        x = x';
    end
    
    if size(grps,1) < size(grps,2)
        grps = grps'; 
    end
    
    %Number of values. 
    n = numel(x);
%% Make figure. 
    %Group numbers. 
    grpNums = unique(grps)';
    jitters = zeros(n,1); 
    
    %Create jitter. 
    c = 1; 
    for g = grpNums
        %Number of elements corresponding to that group. 
        nInGrp = sum(grps==g);
        
        %Jitter!
%         jitters(c:c+nInGrp-1) = g - (sf*randn(nInGrp,1));
        jitters(grps == g) = g - (sf*randn(nInGrp,1));
        
        %Step. 
        c = c+nInGrp; 
    end
    
    %Figure here.
    if isempty(h)
        if isnumeric(position)
            figure('Position',position);
        end
    elseif ~isempty(h)
        axes(h); %Make axes
    end
    hold on;
    hscat = scatter(jitters,x,circleSize,circleColors,'filled');
    alpha(hscat,transparency);
    if plotBox
        hb = boxplot(x,grps,'color',boxColor,'symbol','k','labels',xLabels,...
            'positions',unique(grps));
        boxProps = get(gca,'Children');
        [boxProps(1).Children.LineWidth] = deal(2);
    else 
        hb = nan;
    end
    ylabel(yLabel);
    
    % Plot lines between pairs if specified
    if ~isempty(paired_ind)
        arrayfun(@(a) plot(jitters(paired_ind == a), x(paired_ind == a), ...
            'k-'), unique(paired_ind))
    end
    set(gca,'tickdir','out','box','off');
    
    h = gca;
    
end