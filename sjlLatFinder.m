function latencies = sjlLatFinder(srcRaster,snkRaster)
%
%
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
end