function latencies = sjlLatFinder(srcRaster,snkRaster)
%latencies = sjlLatFinder(srcRaster,snkRaster)
%
%   Finds latencies between transients based on Sam Levy's indexing trick. 
%   Vectorize spike times, cell ID, and lap number and find instances where
%   the diff of cell ID indicates a different neuron and the diff of laps
%   indicates the same lap. 
%
%   INPUTS
%       srcRaster and snkRaster: the first and second rasters (from FT)
%       that you want to analyze.
%
%   OUTPUT
%       latency: latency distribution. 
%

%% Initialize.
   [~,nBins] = size(srcRaster);
    
    %If empty rasters, quit. 
    if ~any(srcRaster(:)) || ~any(snkRaster(:)), return; end;

%% Vectorize spike times then find latencies. 
    [srcLaps,srcTimes] = find(srcRaster);
    [snkLaps,snkTimes] = find(snkRaster); 

    srcTimes = srcTimes + nBins*(srcLaps-1);        %Source spike times.
    snkTimes = snkTimes + nBins*(snkLaps-1);        %Sink spike times.
    
    %Vectorize.
    spktimes = [srcTimes', snkTimes'];
    cellID = [zeros(1,length(srcTimes)), ones(1,length(snkTimes))];
    laps = [srcLaps', snkLaps'];

    %Sort.
    [spktimes,order] = sort(spktimes);
    cellID = cellID(order);
    laps = laps(order);
    
    %Only look at same lap spike timings between distinct cells. Take the
    %diffs of thet vectors then figure out which adjacent values were from
    %the same lap and different neurons. 
    sameLap = diff(laps)==0;           
    diffCells = diff(cellID)==1; 
    latencies = diff(spktimes); 
    latencies = latencies(sameLap & diffCells)./20; 
end