function plotTimeCellGradient2(md,varargin)
%
%
%

%% Parse inputs.
    p = inputParser; 
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('neurons',getTimeCells(md),@(x) isnumeric(x)); 
    p.addParameter('nPoints',50,@(x) isnumeric(x)); 
    p.addParameter('transformation',[]);
    p.addParameter('lineStyle','-',@(x) ischar(x)); 
    p.addParameter('transparent',false,@(x) islogical(x));
    
    p.parse(md,varargin{:});
    neurons = p.Results.neurons; 
    if size(neurons,1) > size(neurons,2), neurons = neurons'; end
    nPoints = p.Results.nPoints;
    transformation = p.Results.transformation; 
    lineStyle = p.Results.lineStyle; 
    transparent = p.Results.transparent; 
    
%%
    cd(md.Location); 
    
    %Get time of peak. 
    [~,peak] = getTimePeak(md); 
    
    %Time vector. 
    t = linspace(0,10,nPoints); 
    c = colormap(jet(nPoints)); 
    
    %PlotNeurons(md,1:nNeurons,[.7 .7 .7],1); 
    nNeurons = length(peak);
    colors = zeros(nNeurons,3);
    for n=neurons
        closestTime = findclosest(peak(n),t); 
        colors(n,:) = c(closestTime,:); 
    end
            
    hold on;
    PlotNeurons(md,neurons,colors,2,'transformation',transformation,...
        'lineStyle',lineStyle,'transparent',transparent); 
    load('MovieDims.mat','Xdim','Ydim');
    xlim([0 Ydim]); ylim([0 Xdim]); 
end