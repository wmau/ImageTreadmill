function NeighborhoodTimeDistance(md,thresholds)
%NeighborhoodTimeDistance(md)
%
%   Investigate whether time cells cluster based on their time interval
%   coding proximity. Method: Define a distance threshold. For each time
%   cell, look for other time cells in the vicinity of that threshold and
%   get their average spike times. 

%%
    cd(md.Location); 
    load('TimeCells.mat','TimeCells');
    load('TemporalInfo.mat','sig');
    TimeCells = intersect(TimeCells,find(sig)); 
    nTCs = length(TimeCells);
    nDistances = length(thresholds);
    B = 1000;
    
%% Get centroids. 
    centroids = getNeuronCentroids(md,'neurons',TimeCells);
    centroids = centroids(~isnan(centroids(:,2)),:);
    
%% Get anatomical distances.
    D = nan(nTCs);
    for n1=1:nTCs
        %Centroid for cell 1.
        x1 = centroids(n1,1);
        y1 = centroids(n1,2); 
        
        for n2=n1+1:nTCs
            %centroid for cell 2.
            x2 = centroids(n2,1);
            y2 = centroids(n2,2);

            %Anatomical distance. 
            D(n1,n2) = sqrt((x2-x1)^2 + (y2-y1)^2);
        end
    end
    
%% Get time peaks. 
    [~,t] = getTimePeak(md);
    t = t(TimeCells);
%     [~,t] = PastalkovaPlot(md,false);
%     t = t./max(t);
    
%% 
    %tEdges = 0:.05:1;
    tEdges = 0:.5:10;
    p = zeros(length(tEdges),nDistances);
    [neighbors,dT] = deal(cell(nTCs,nDistances));
    avg_dT = zeros(1,nDistances);
    for d=1:nDistances
        thisThresh = thresholds(d);
        
        for c=1:nTCs
            neighbors{c,d} = find(D(c,c+1:end) <= thisThresh);
 
            dT{c,d} = abs(t(neighbors{c,d}) - t(c));
        end
        
        dTsforthisThresh = cell2mat(dT(:,d));
        
        avg_dT(d) = mean(dTsforthisThresh);
        %p(:,d) = histc(dTsforthisThresh,tEdges)./length(dTsforthisThresh);
    end
    
%     figure; 
%     surface(thresholds,tEdges,p);
%     view(3);
%     shading interp;

%% Shuffle neuron identities in time peak then rerun analysis. 
    [rNeighbors,rdT] = deal(cell(nTCs,nDistances));
    rAvg_dT = zeros(B,nDistances);
    prog = ProgressBar(B);
    for i=1:B
        %Shuffle time peaks.
        rT = t(randperm(length(t)));
        
        for d=1:nDistances
            thisThresh = thresholds(d);

            for c=1:nTCs
                rNeighbors{c,d} = find(D(c,c+1:end) <= thisThresh);

                rdT{c,d} = abs(rT(rNeighbors{c,d}) - t(c));
            end

            dTsforthisThresh = cell2mat(rdT(:,d));

            rAvg_dT(i,d) = mean(dTsforthisThresh);
        end
        
        prog.progress;
    end
    prog.stop;
    
    rAvg_dT = sort(rAvg_dT);
    surrogate = mean(rAvg_dT); 
    ci(:,1) = surrogate - rAvg_dT(round(.01*B),:);
    ci(:,2) = -(surrogate - rAvg_dT(round(.99*B),:));
    figure;
    plot(thresholds,avg_dT,'k','linewidth',3);
    hold on;
    l = boundedline(thresholds,surrogate,ci,'alpha');
    l.Color = [.5 .5 1]; l.LineStyle = '--';
    xlabel('Distance Threshold [\mum]');
    ylabel('Rank Distance');
    axis tight; 
    set(gca,'tickdir','out');
end