function [d,td] = TimeTopology(md)
%
%
%

%%
    cd(md.Location);
    load('TimeCells.mat','TimeCells');
    
    peaks = getTimePeak(md);
    centroids = getNeuronCentroids(md,'neurons',TimeCells);
    nTCs = length(TimeCells);
    nNeurons = length(peaks);
    
    d = nan(nNeurons);
    td = nan(nNeurons);
    
    for i=TimeCells'
        x1 = centroids(i,1);
        y1 = centroids(i,2); 
        
        pool = TimeCells(TimeCells > i);
        for j=pool'
            x2 = centroids(j,1); 
            y2 = centroids(j,2);

            d(i,j) = sqrt((x2-x1)^2 + (y2-y1)^2);
            td(i,j) = peaks(i) - peaks(j);

        end
    end
    
end