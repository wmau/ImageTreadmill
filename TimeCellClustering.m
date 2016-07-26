function p= TimeCellClustering(md)
%
%
%

%%
    B = 1000;
    cd(md.Location);
    
    load('TimeCells.mat','TimeCells');
    nTCs = length(TimeCells);
    peaks = getTimePeak(md);
    peaks = peaks(TimeCells);
    
    centroids = getNeuronCentroids(md);
    centroids = centroids(TimeCells,:);
    
    [~,order] = sort(peaks);
    
    peaks = peaks(order); 
    centroids = centroids(order,:);
    
    dists = nan(1,nTCs-1);
    null = nan(B,nTCs-1);
    nullmeans = nan(1,B);
    for j=1:nTCs-1
        dists(j) = pdist(centroids(j:j+1,:),'euclidean');
    end
    distmean = median(dists);
    
    p = ProgressBar(B);
    for i=1:B
        shuffled = centroids(randperm(nTCs),:);
        for j=1:nTCs-1
            null(i,j) = pdist(shuffled(j:j+1,:),'euclidean');     
        end
        nullmeans(i) = median(null(i,:));
        
        p.progress;
    end
    p.stop;
    
    p = sum(distmean > nullmeans)/B;
    
    figure;
    histogram(nullmeans,'normalization','probability');
    yl = ylim;
    hold on;
    line([distmean distmean],[0 yl(2)],'color','r');
    hold off
    set(gca,'ticklength',[0 0]);
    xlabel('Mean Distances [pixels]');
    ylabel('Proportion of Iterations');
    
    figure;
    histogram(dists,30,'facecolor','k');
    xlabel('Distances Between Adjacent Time Cells');
    ylabel('Count');
    set(gca,'ticklength',[0 0]);
    
end