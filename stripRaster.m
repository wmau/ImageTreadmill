function [lagRaster,closest] = stripRaster(lagRaster,leadRaster) 
%lagRaster = stripRaster(lagRaster,leadRaster)
%
%   Takes lagRaster (raster of neuron 1, preceding neuron 2) and takes out
%   all spikes except for those immediately preceding a spike in
%   leadRaster. 
%
%   INPUTS
%       lagRaster: logical (trial x timestamp) array denoting spike times
%       for neuron 1. 
%
%       leadRaster: same as lagRaster, but for neuron 2. 
%
%   OUTPUT
%       lagRaster: new lagRaster, but with only those responses immediately
%       preceding neuron 2. 
%

%% Main body. 
    %Preallocate. 
    temp = false(size(lagRaster)); 
    closest = [];
    
    %Find indices of spikes. 
    [lagLap,lagT] = find(lagRaster); 
    [leadLap,leadT] = find(leadRaster); 
    
    %For each lap where both neurons were active...
    for l=intersect(leadLap,lagLap)'
        %Get all spike times of neuron 2 on that lap. 
        leadOnset = leadT(leadLap==l);
        
        %For each spike of neuron 2 on lap l...
        for o=leadOnset'
            %Get all spike times of neuron 1 on that lap. 
            lagOnsets = lagT(lagLap==l);
            
            %Filter out those that occurred after neuron 2.
            lagOnsets = lagOnsets(lagOnsets<=o);
      
            %Get closest spike. 
            bingo = lagOnsets(findclosest(o,lagOnsets));
            closest = [closest; bingo - o];
        
            %Tick. 
            temp(l,bingo) = true; 
        end
    end
    
    %Divide by frame rate. 
    closest = closest./20; 
    
    %Replace lagRaster with processed 
    lagRaster = temp; 
end