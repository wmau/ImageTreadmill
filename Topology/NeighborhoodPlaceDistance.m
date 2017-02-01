function NeighborhoodPlaceDistance(md,thresholds)
%
%
%

%%
    cd(md.Location); 
    load(fullfile(pwd,'Placefields.mat'),'pval');
    load(fullfile(pwd,'PlacefieldStats.mat'),'PFnHits','bestPF');
    load(fullfile(pwd,'SpatialInfo.mat'),'MI'); 
    idx = sub2ind(size(PFnHits),1:size(PFnHits,1),bestPF');
    PCcrit = .01;
    PCs = find(pval<PCcrit & MI'>0 & PFnHits(idx)>4); 
    nPCs = length(PCs);
    nDistances = length(thresholds);
    B = 1000;
    
%% 
    centroids = getNeuronCentroids(md,'neurons',PCs);
    centroids = centroids(~isnan(centroids(:,2)),:);
    
    D = nan(nPCs);
    for n1=1:nPCs
        %Centroid for cell 1.
        x1 = centroids(n1,1);
        y1 = centroids(n1,2); 
        
        for n2=n1+1:nPCs
            %centroid for cell 2.
            x2 = centroids(n2,1);
            y2 = centroids(n2,2);

            %Anatomical distance. 
            D(n1,n2) = sqrt((x2-x1)^2 + (y2-y1)^2);
        end
    end
    
%% 
    [~,~,~,order] = LinearizedPFs_treadmill(md);
    order = order./max(order); 
    
%% 
    pEdges = 0:.05:1;
    p = zeros(length(pEdges),nDistances);
    [neighbors,dX] = deal(cell(nPCs,nDistances));
    avg_dX = zeros(1,nDistances);
    for d=1:nDistances
        thisThresh = thresholds(d);
        
        for c=1:nPCs
            neighbors{c,d} = find(D(c,c+1:end) <= thisThresh);
 
            dX{c,d} = abs(order(neighbors{c,d}) - order(c));
        end
        
        dXsforthisThresh = cell2mat(dX(:,d));
        
        avg_dX(d) = mean(dXsforthisThresh);
        p(:,d) = histc(dXsforthisThresh,pEdges)./length(dXsforthisThresh);
    end
    
%% 
    [rNeighbors,rdX] = deal(cell(nPCs,nDistances));
    rAvg_dX = zeros(B,nDistances);
    prog = ProgressBar(B);
    for i=1:B
        %Shuffle time peaks.
        rX = order(randperm(length(order)));
        
        for d=1:nDistances
            thisThresh = thresholds(d);

            for c=1:nPCs
                rNeighbors{c,d} = find(D(c,c+1:end) <= thisThresh);

                rdX{c,d} = abs(rX(rNeighbors{c,d}) - order(c));
            end

            dXsforthisThresh = cell2mat(rdX(:,d));

            rAvg_dX(i,d) = mean(dXsforthisThresh);
        end
        
        prog.progress;
    end
    prog.stop;
    
    rAvg_dX = sort(rAvg_dX);
    surrogate = mean(rAvg_dX); 
    ci(:,1) = surrogate - rAvg_dX(round(.01*B),:);
    ci(:,2) = -(surrogate - rAvg_dX(round(.99*B),:));
    figure;
    plot(thresholds,avg_dX,'k','linewidth',3);
    hold on;
    l = boundedline(thresholds,surrogate,ci,'alpha');
    l.Color = [.5 .5 1]; l.LineStyle = '--';
    xlabel('Distance Threshold [\mum]');
    ylabel('Rank Distance');
    axis tight; 
    set(gca,'tickdir','out');
end