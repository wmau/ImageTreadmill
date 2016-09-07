function [p,el] = trialShuffleKSTest(graphData,neuron,varargin)
%[p,el] = trialShuffleKSTest(graphData,neuron,varargin)
%
%   Perform randomized trial shuffles for treadmill run epochs where two
%   specified cells are active. Then compare their latency distribution to
%   the latency distribution after the trial shuffle. 
%
%   INPUTS
%       graphData: Output from MakeGraphv4.
%
%       neuron: Putative target neuron.
%
%       varargin (optional)
%           edgelist: Vector of putative trigger neurons. If not provided,
%           will use adjacency matrix in graphData.
%
%   OUTPUTS
%       p: 1xE (E = number of putative triggers) vector with p-values for
%       KS tests.
%
%       el: 1xE vector, edge list (same as input if specified). 
%

%% Set up. 
    p = inputParser;
    p.addRequired('graphData',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x)); 
    p.addParameter('edgelist',find(graphData.A(:,neuron))',@(x) isnumeric(x));
    p.parse(graphData,neuron,varargin{:});
    
    %Assign values. 
    el = p.Results.edgelist; 
    neuron = p.Results.neuron;
    
    %Number of neurons that presumably connect to neuron. 
    nTriggers = length(el);
    
    %Number of trial shuffles. 
    B=1000;
    
    %Get MD entry. 
    md = findMDfromGraphData(graphData); 
  
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
    
    %Build raster for target neuron. 
    targetRaster = buildRaster(inds,FT,neuron);
    
    %Preallocate things. 
    c=1;                        %Counter.
    p = ones(1,nTriggers);      %p-values.
    
%% Do trial shuffle. 
    %Get treadmill spread for each connection.
    [~,~,treadmillSpread,TMAlignedOnsets] = SpreadRatio(md,graphData,neuron,'inds',inds);
    
    for e=el 
        %If treadmill-target latency is variable...
        if treadmillSpread(c) > 0 && mode(TMAlignedOnsets{c})~=0
            %Build the tick raster for the trigger neuron. 
            triggerRaster = buildRaster(inds,FT,e); 
            
            %Get trigger-target latencies. 
            [~,latencies] = stripRaster(triggerRaster,targetRaster); 
            null = [];
            
            %Find intersect of set of laps. 
            [targetLaps,~] = find(targetRaster); 
            [triggerLaps,~] = find(triggerRaster); 
            laps = intersect(targetLaps,triggerLaps); 
            nLaps = length(laps); 
            
            %Only look at those laps. 
            tempTarget = targetRaster(laps,:);
            tempTrigger = triggerRaster(laps,:); 
            
            %Shuffle trials then append latencies to null distribution.
            for i=1:B
                shuffledTrigger = tempTrigger(randperm(nLaps),:);
                [~,shuff_latencies] = stripRaster(shuffledTrigger,tempTarget);  
                null = [null; shuff_latencies]; 
            end
            
            %Do KS test after B shuffles. 
            [~,p(c)] = kstest2(latencies,null); 
            
            c=c+1; 
        end
    end
    
    %Round to get rid of round-off error. 
    p = round(p,4);
end