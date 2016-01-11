function clusterPastalkova(animal,date,session,T,vv,clusters)
%
%
%

%% 
    ChangeDirectory(animal,date,session); 
    try
        load('TimeCells.mat');
    catch
        [TimeCells,ratebylap,curves,delays,x,y,time_interp] = FindTimeCells(animal,date,session,T); 
    end
    
    %Useful variables. 
    [nLaps,nBins,nNeurons] = size(ratebylap); 
    nClusters = length(clusters); 
    
    figure;
    for i=1:nClusters       
        %Find neurons that live in this cluster. 
        thisCluster = clusters(i);
        inCluster = vv==thisCluster;
        tilemat = cell2mat(curves.tuning(inCluster)); 
        
        %Find the peak, normalize, and sort. 
        [peaks,peakInds] = max(tilemat,[],2);
        normtilemat = tilemat./repmat(peaks,1,nBins);
        [~,order] = sort(peakInds);
        sortedMat = normtilemat(order,:); 
        
        %Plot.
        subplot(nClusters,1,i); 
        imagesc([0:T],[1:5:sum(inCluster)],sortedMat);
            colormap gray; ylabel('Neurons');
            if i==nClusters
                xlabel('Time [s]'); 
            end
            
    end
    
end