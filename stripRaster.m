function [trigRaster,closest,lap] = stripRaster(trigRaster,targRaster) 
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
    temp = false(size(trigRaster)); 
    closest = [];
    lap = [];
    
    %Find indices of spikes. 
    [trigLap,trigT] = find(trigRaster); 
    [targLap,targT] = find(targRaster); 
    
    %For each lap where both neurons were active...
    for l=intersect(targLap,trigLap)'
        %Get all spike times of neuron 2 on that lap. 
        targOnset = targT(targLap==l);
       
        %For each spike of neuron 2 on lap l...
        for o=targOnset'
            %if l==10, keyboard; end;
            %Get all spike times of neuron 1 on that lap. 
            trigOnsets = trigT(trigLap==l);
            
            %Filter out those that occurred after neuron 2.
            trigOnsets = trigOnsets(trigOnsets<=o);
      
            %Get closest spike. 
            bingo = trigOnsets(findclosest(o,trigOnsets));

            if ~isempty(bingo)
                closest = [closest; bingo - o];
                ind = find(trigT <= bingo & trigLap==l);
                trigT(ind) = []; 
                trigLap(ind) = []; 
                
                lap = [lap; l]; 
            end;
        
            %Tick. 
            temp(l,bingo) = true; 
        end
    end
    
    %Divide by frame rate. 
    closest = closest./20; 
    
    %Replace lagRaster with processed 
    trigRaster = temp; 
end