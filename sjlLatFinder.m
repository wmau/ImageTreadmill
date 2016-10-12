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
    latencies = [];
    nLaps = size(srcRaster,1);
    
    %If empty rasters, quit. 
    if ~any(srcRaster(:)) || ~any(snkRaster(:)), return; end;

%% Vectorize spike times then find latencies. 
    spktimes = [];
    cellID = [];
    laps = [];
    for l=1:nLaps
        srcSpks = find(srcRaster(l,:));     %Source spike times.
        snkSpks = find(snkRaster(l,:));     %Sink spike times.
        
        %Vectorize spike times, cell ID, and lap number.
        spktimes = [spktimes, srcSpks, snkSpks];
        cellID = [cellID, zeros(1,length(srcSpks)), ones(1,length(snkSpks))];   
        laps = [laps, l*ones(1,length(srcSpks)), l*ones(1,length(snkSpks))];
    end
    
    %Only look at same lap spike timings between distinct cells. Take the
    %diffs of thet vectors then figure out which adjacent values were from
    %the same lap and different neurons. 
    sameLap = diff(laps)==0;           
    diffCells = diff(cellID)==1; 
    latencies = diff(spktimes); 
    latencies = latencies(sameLap & diffCells)./20; 
    latencies = latencies(latencies > 0);
end