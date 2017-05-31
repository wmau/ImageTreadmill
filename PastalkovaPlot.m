function [sortedPastalkova,order,peakInds,tilemat] = PastalkovaPlot(MD,varargin)
%sortedPastalKova = PastalkovaPlot(animal,date,session,T)
%
%   Makes a plot that shows response curves that tile the delay. 

%% Make a plot that shows response curves that tile the delay. 
    p = inputParser;
    p.addRequired('MD',@(x) isstruct(MD));
    p.addParameter('plotit',true,@(x) islogical(x)); 
    p.addParameter('TimeCells',getTimeCells(MD),@(x) isnumeric(x));
    p.addParameter('order',false);
    
    p.parse(MD,varargin{:});
    
    TimeCells = p.Results.TimeCells;
    plotit = p.Results.plotit;
    order = p.Results.order;
    
    cd(MD.Location); 
    load('TimeCells.mat','curves','T');
      
    %Concatenate all time cell tuning curves.
    tilemat = cell2mat(curves.tuning(TimeCells));
    nBins = size(tilemat,2); 
    
    %Find the peaks.
    [peaks,peakInds] = max(tilemat,[],2);
    if ~order, [peakInds,order] = sort(peakInds); 
    else, peakInds = peakInds(order); end
    
    %Normalize.
    normtilemat = tilemat./repmat(peaks,1,nBins);
    
    %Sort.
    tilemat = tilemat(order,:);
    sortedPastalkova = normtilemat(order,:);
    
    %Plot. 
    if plotit
    imagesc([0:T],[1:length(TimeCells)],sortedPastalkova);
        colormap gray;
        set(gca,'ydir','reverse','ytick',[1,length(TimeCells)]); axis tight;
        xlabel('Time [s]'); ylabel('Neurons'); 
    end

end