function [p,el] = testSpreadRatio(md,graphData,neuron,varargin)
%[p,el] = testSpreadRatio(md,graphData,neuron,varargin)
%
%   Statistical test to determine whether or not the "spread ratio" of a
%   cell pair is significant. The spread ratio is the median absolute
%   deviation of neuron two's temporal relationship to neuron one divided
%   by neuron two's temporal relationship with treadmill onset (only laps
%   where both neurons were active). Performs B trial identity shuffles and
%   checks whether the temporal relationship between the cell pairs are
%   still consistent. 
%
%   INPUTS
%       md: Session entry
%
%       graphData: Output from MakeGraphvX. 
%
%       neuron: Neuron two. 
%
%       varargin valid inputs...
%           edgelist: Vector containing neurons whose connections to neuron
%           two you want to examine. 
%
%   OUTPUTS
%       p: Vector containing p-values for each cell pair temporal
%       consistency. 
%
%       el: Edge list vector. Cells that presumably connect to the "neuron"
%       argument. 
%

%% Grab inputs.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(graphData.A(:,neuron))',@(x) isnumeric(x));
    p.parse(md,graphData,neuron,varargin{:});
    
    el = p.Results.edgelist; 
    md = p.Results.md; 
    neuron = p.Results.neuron;
    
    %Number of neurons that presumably connect to neuron. 
    nInitiators = length(el);
    
    B=1000;
%% Set up.
    %Change directory and load initial variables. 
    cd(md.Location); 
    load('TimeCells.mat','T','TodayTreadmillLog'); 
    load('Pos_align.mat','FT');
    load('ProcOut.mat','NumNeurons');
    complete = TodayTreadmillLog.complete;
    inds = TodayTreadmillLog.inds; 
    
    %Get treadmill run indices. 
    inds = inds(find(complete),:);  %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1; %Consistent length.   
    nLaps = sum(complete); 
    
    %Build raster for neuron 2. 
    leadRaster = buildRaster(inds,FT,neuron);
    
    %Preallocate things. 
    c=1; 
    ratio = zeros(1,nInitiators);
    null = zeros(1,B);
    p = zeros(1,nInitiators);
    
%% For each neuron in the edge list, test consistency by doing trial shuffles. 
    for e=el
        %Build the tick raster for neuron 1. 
        lagRaster = buildRaster(inds,FT,e);        
        [immediateRaster,d] = stripRaster(lagRaster,leadRaster);

        %Get spread. 
        cellOffsetSpread = mad(d,1);
        
        %Only look at laps where both neurons were active. 
%         bothActiveLaps = find(any(immediateRaster,2)); 
%         TMalignedOnsets = [];
%         for l=bothActiveLaps'
%             %Get the onset times of each neuron. 
%             TMalignedOnsets = [TMalignedOnsets find(leadRaster(l,:))];     
%         end
        [~,TMalignedOnsets] = find(leadRaster);

        %Divide by frame rate . 
        TMalignedOnsets = TMalignedOnsets./20;

        %Spread of responses relative to treadmill start. 
        treadmillOffsetSpread = mad(TMalignedOnsets,1);

        %Ratio between cell-to-cell vs cell-to-treadmill.
        ratio(c) = cellOffsetSpread / treadmillOffsetSpread;
        
        %Shuffle trial identity and find the new spread ratio. 
        for i=1:B
            shuffled = lagRaster(randperm(nLaps),:);
            [~,dB] = stripRaster(shuffled,leadRaster);
            dBSpread = mad(dB,1);
            null(i) = dBSpread./treadmillOffsetSpread;        
        end
        
        %P-value, number of times where the real temporal relationship
        %spread was wider than the shuffled. 
        p(c) = sum(ratio(c) > null) / B;
      
        c=c+1;  %counter.
    end
    
end