function boxScatterplot(x,grps,varargin)
%
%
%

%%
    p = inputParser; 
    p.addRequired('x',@(x) isnumeric(x));
    p.addRequired('grps',@(x) isnumeric(x));
    p.addParameter('xLabels',@(x) ischar(x));
    p.addParameter('yLabel',@(x) ischar(x)); 
    p.addParameter('boxColor','k',@(x) ischar(x) || isnumeric(x));
    p.addParameter('circleSize',10,@(x) isnumeric(x)); 
    p.addParameter('circleColors',[.7 .7 .7],@(x) ischar(x) || isnumeric(x));
    p.addParameter('transparency',.5,@(x) isscalar(x)); 
    p.addParameter('sf',.05,@(x) isscalar(x));
    p.addParameter('position',[520 378 560 420],@(x) isnumeric(x)); 
    
    p.parse(x,grps,varargin{:});
    xLabels = p.Results.xLabels; 
    yLabel = p.Results.yLabel; 
    boxColor = p.Results.boxColor;
    circleSize = p.Results.circleSize;
    transparency = p.Results.transparency;
    circleColors = p.Results.circleColors;
    sf = p.Results.sf; 
    position = p.Results.position;
    
    %Turn into column.
    if size(x,1) < size(x,2) 
        x = x';
    end
    
    if size(grps,1) < size(grps,2)
        grps = grps'; 
    end
    
    n = numel(x);
%% 
    grpNums = unique(grps)';
    jitters = zeros(n,1); 
    
    i = 1; 
    c = 1; 
    for g = grpNums
        nInGrp = sum(grps==g);
        jitters(c:c+nInGrp-1) = i - (sf*randn(nInGrp,1));
        
        i = i+1;
        c = c+nInGrp; 
    end
    
    figure('Position',position); hold on;
    scat = scatter(jitters,x,circleSize,circleColors,'filled');
    alpha(scat,transparency);
    boxplot(x,grps,'color',boxColor,'symbol','k','labels',xLabels);
    ylabel(yLabel);
    set(gca,'ticklength',[0 0]);
    
end