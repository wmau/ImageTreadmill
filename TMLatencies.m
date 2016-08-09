function TMAlignedOnsets = TMLatencies(immRaster,targRaster)
%TMAlignedOnsets = TMLatencies(immRaster,targRaster)
%
%   Gets the treadmill-target latencies for laps where both the trigger and
%   target cells are active. 
%
%   INPUTS
%       immRaster: trigger raster with all spikes removed except for those
%       immediately preceding spikes from the target. 
%
%       targRaster: target cell raster.
%
%   OUTPUT
%       TMAlignedOnsets: Treadmill-target latencies (in seconds). 
%

%% Main body. 
    bothActiveLaps = find(any(immRaster,2)); 
    TMAlignedOnsets = [];
    
    %Only on laps where both are active. 
    for l=bothActiveLaps'
        %Get the onset times of each neuron. 
        TMAlignedOnsets = [TMAlignedOnsets find(targRaster(l,:))];    
    end
    
    %Divide by frame rate then subtract by 1 frame. 
    TMAlignedOnsets = TMAlignedOnsets ./ 20 - 0.05; 
    
end