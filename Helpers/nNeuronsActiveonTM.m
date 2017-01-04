function [n,neurons] = nNeuronsActiveonTM(md)
%[n,neurons] = nNeuronsActiveonTM(md)
%
%   Simply gets the number of neurons that are active on the treadmill for
%   more laps than the critical number (a quarter of the total laps). 
%
%   INPUT
%       md: session entry.
%
%   OUTPUTS
%       n: number of neurons. 
%
%       neurons: neurons active. 
%

%% Set up.
    cd(md.Location);
    load('TimeCells.mat','TodayTreadmillLog','T');
    
    %Get aligned FTs. 
    try
        load('Pos_align.mat','PSAbool'); 
    catch
        load('FinalOutput.mat','PSAbool');
        [~,~,~,PSAbool] = AlignImagingToTracking(md.Pix2CM,PSAbool,0); 
    end
    
    nNeurons = size(PSAbool,1);
    inds = TodayTreadmillLog.inds;
    inds = inds(find(TodayTreadmillLog.complete),:);        %Only completed runs. 
    inds(:,2) = inds(:,1) + 20*T-1;                         %Consistent length.
    nRuns = sum(TodayTreadmillLog.complete); 
    critLaps = 0.25*nRuns; 
    
    raster = cell(1,nNeurons);
    %Build all the rasters.
    for n=1:nNeurons
        raster{n} = buildRaster(inds,PSAbool,n);    
    end
    
    %Find neurons active for more than the required number of laps. 
    nLapsActive = cell2mat(cellfun(@(x) sum(any(x,2)),raster,'unif',0));
    neurons = find(nLapsActive > critLaps);
    n = length(neurons);
   
end