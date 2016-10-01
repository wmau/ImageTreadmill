function [ratio,cellSpread,treadmillSpread,TMAlignedOnsets] = SpreadRatio(md,A,neuron,varargin)
%[ratio,cellSpread,treadmillSpread,TMAlignedOnsets] = SpreadRatio(md,graphData,neuron,varargin)
%
%

%% Grab inputs.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('A',@(x) isnumeric(x) || islogical(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(A(:,neuron))',@(x) isnumeric(x));
    p.addParameter('inds',false);
    p.parse(md,A,neuron,varargin{:});
    
    el = p.Results.edgelist; 
    md = p.Results.md; 
    neuron = p.Results.neuron;
    inds = p.Results.inds;
    
    %Number of neurons that presumably connect to neuron. 
    nInitiators = length(el);
    
%%
    cd(md.Location); 
    load('Pos_align.mat','FT');
    
    if ~isnumeric(inds)
        load('TimeCells.mat','T','TodayTreadmillLog');
        
        complete = TodayTreadmillLog.complete;
        inds = TodayTreadmillLog.inds; 

        %Get treadmill run indices. 
        inds = inds(find(complete),:);  %Only completed runs. 
        inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    end
    
    %Build raster for neuron 2. 
    targetRaster = buildRaster(inds,FT,neuron);

    %Treadmill-aligned onsets. 
    %[~,TMAlignedOnsets] = find(leadRaster); 

    %Divide by frame rate.
    %TMAlignedOnsets = TMAlignedOnsets./20;

    %Get cell-treadmill latency spread. 
    %treadmillSpread = mad(TMAlignedOnsets,1); 
%%
    c=1;
    mm=1;
    ratio = nan(1,nInitiators); 
    cellSpread = nan(1,nInitiators); 
    TMAlignedOnsets = cell(1,nInitiators);
    treadmillSpread = nan(1,nInitiators);
    for e=el
        %Build the tick raster for neuron 1. 
        triggerRaster = buildRaster(inds,FT,e);        
        [immediateRaster,latencies] = stripRaster(triggerRaster,targetRaster);

        %Get treadmill-target latencies.
        TMAlignedOnsets{c} = TMLatencies(immediateRaster,targetRaster);
        
        %Get treadmill latency spread. 
        treadmillSpread(c) = mad(TMAlignedOnsets{c},mm); 

        %Get cell pairing latency spread. 
        cellSpread(c) = mad(latencies,mm);

        %Get cell pairing to cell-treadmill latency spread ratio.
        ratio(c) = cellSpread(c) / treadmillSpread(c); 

        c = c+1; 
    end
    
end