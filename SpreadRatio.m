function [ratio,cellSpread,treadmillSpread,TMAlignedOnsets] = SpreadRatio(md,graphData,neuron,varargin)
%[ratio,cellSpread,treadmillSpread,TMAlignedOnsets] = SpreadRatio(md,graphData,neuron,varargin)
%
%

%% Grab inputs.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(graphData.A(:,neuron))',@(x) isnumeric(x));
    p.addParameter('inds',false);
    p.parse(md,graphData,neuron,varargin{:});
    
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
    leadRaster = buildRaster(inds,FT,neuron);

    %Treadmill-aligned onsets. 
    %[~,TMAlignedOnsets] = find(leadRaster); 

    %Divide by frame rate.
    %TMAlignedOnsets = TMAlignedOnsets./20;

    %Get cell-treadmill latency spread. 
    %treadmillSpread = mad(TMAlignedOnsets,1); 
%%
    c=1;
    ratio = nan(1,nInitiators); 
    cellSpread = nan(1,nInitiators); 
    TMAlignedOnsets = cell(1,nInitiators);
    treadmillSpread = nan(1,nInitiators);
    for e=el
        %Build the tick raster for neuron 1. 
        lagRaster = buildRaster(inds,FT,e);        
        [immediateRaster,d] = stripRaster(lagRaster,leadRaster);

        %Only look at laps where both neurons were active. 
        bothActiveLaps = find(any(immediateRaster,2)); 

        %Only look at laps where both neurons were active. 
        for l=bothActiveLaps'
            %Get the onset times of each neuron. 
            TMAlignedOnsets{c} = [TMAlignedOnsets{c} find(leadRaster(l,:))];     
        end
        TMAlignedOnsets{c} = TMAlignedOnsets{c}./20;

        %Get treadmill latency spread. 
        treadmillSpread(c) = mad(TMAlignedOnsets{c},1); 

        %Get cell pairing latency spread. 
        cellSpread(c) = mad(d,1);

        %Get cell pairing to cell-treadmill latency spread ratio.
        ratio(c) = cellSpread(c) / treadmillSpread(c); 

        c = c+1; 
    end
    
end