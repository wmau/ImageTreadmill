function [d,td] = TimeCellTopology(md)
%[d,td] = TimeCellTopology(md)
%
%   One version of time cell clustering analysis. Tries to answer the
%   question "Are time cells that code a similar interval closer to each
%   other?" Method: rank time cells based on their tuning curve peaks and
%   get all pairwise rank differences and anatomical distances. Then
%   correlate. 
%

%%
    cd(md.Location);
    load('TimeCells.mat','TimeCells','T'); 
    load('TemporalInfo.mat','sig');
    TimeCells = intersect(find(sig),TimeCells);
    nTCs = length(TimeCells);
    
    [~,order] = PastalkovaPlot(md,T,false);
    centroids = getNeuronCentroids(md,'neurons',TimeCells);
    
    [d,td] = deal(nan(nTCs));
    for i=1:nTCs
        %Centroid for cell 1.
        x1 = centroids(TimeCells(i),1);
        y1 = centroids(TimeCells(i),2); 
        
        for j=i+1:nTCs
            %Centroid for cell 2.
            x2 = centroids(TimeCells(j),1); 
            y2 = centroids(TimeCells(j),2);

            %Anatomical distance. 
            d(i,j) = sqrt((x2-x1)^2 + (y2-y1)^2);
            
            %Sequence rank difference. 
            td(i,j) = order(i) - order(j);

        end
    end
    
    d = d(~isnan(d));
    td = abs(td(~isnan(td)));
    
    dots = scatter(td,d,10,'filled');
    alpha(dots,.2); 
    lsline;
    xlabel('Rank distance'); 
    ylabel('Anatomical distance [\mum]');
    
end