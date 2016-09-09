function raster = buildRaster(inds,FT,neuron) 
%raster = buildRaster(nTrials,inds,FT,neuron)
%
%   Builds raster of calcium transient onset times at the native sampling
%   rate of 20 Hz.
%
%   INPUTS
%       inds: trialsx2 matrix of column indices of FT specifying the
%       timestamps at which to build rasters.
%
%       FT: From T2output.
%
%       neuron: Scalar, row index of FT specifying which neuron.
%

%% Build raster. 
    nTrials = size(inds,1);                             %Number of trials.
    treadmilldurations = unique(diff(inds,[],2)); 
    if length(treadmilldurations) > 1
        disp('Multiple treadmill durations found:'); 
        disp(num2str(treadmilldurations)); 
        disp('Equivalent to:'); 
        disp([num2str([treadmilldurations./20]'), ' seconds.']);
        disp('Correcting using 200 frames (10 seconds)'); 
        
        inds(:,2) = inds(:,1) + 199;
        raster = zeros(nTrials,200);
    else
        raster = zeros(nTrials,treadmilldurations+1);  %Preallocate.
    end
    for t=1:nTrials
        lapRaster = [0 FT(neuron,inds(t,1):inds(t,2))]; %Take column difference of binary FT. 
        raster(t,:) = diff(lapRaster);          
    end
    raster(raster == -1) = 0;                           %Erase offsets. 
    raster = logical(raster);                           %Turn into logical. 
end