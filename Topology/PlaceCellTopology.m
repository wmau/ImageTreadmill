function PlaceCellTopology(md)
%
%
%

%%
    cd(md.Location);
    PCs = getPlaceCells(md,.01);
    nPCs = length(PCs);
    
    
    [~,~,~,order] = LinearizedPFs_treadmill(md,'plotit',false);
    order = order./max(order);
    
    centroids = getNeuronCentroids(md,'neurons',PCs);
    
    [d,pd] = deal(nan(nPCs));
    for i=1:nPCs
        %Centroid for cell 1.
        x1 = centroids(PCs(i),1);
        y1 = centroids(PCs(i),2);
        
        for j=i+1:nPCs
            %Centroid for cell 2.
            x2 = centroids(PCs(j),1);
            y2 = centroids(PCs(j),2);
            
            %Anatomical distance.
            d(i,j) = sqrt((x2-x1)^2 + (y2-y1)^2);
            
            %Sequence rank difference.
            try
            pd(i,j) = order(i) - order(j);
            catch, keyboard; end 
            
        end
        
    end
    
    d = d(~isnan(d));
    pd = abs(pd(~isnan(pd)));
    [R,p] = corr(pd,d);
    
    dots = scatter(pd,d,10,'filled');
    alpha(dots,.2); 
    lsline;
    xlabel('Rank distance'); 
    ylabel('Anatomical distance [microns]');
    title(['p = ',num2str(p)]);
end
    