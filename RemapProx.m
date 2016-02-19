function RemapProx(tday,yday,MAPlocation,Ts)
%
%
%

%% 
    %Make time vector. 
    tResolution = 0.25; 
    nBins = Ts(2)/tResolution; 
    ynBins = Ts(1)/tResolution;
    t = linspace(0,Ts(2),nBins)';
    yt = linspace(0,Ts(1),ynBins)'; 
    
    %Find the neurons in tday that are either new time cells or have
    %remapped since yday. 
    [remapped,MAP,MAPcols,TIMECELLS,CURVES] = whichTimeCellsAreNew(tday,yday,MAPlocation,Ts);
    
    %Get the indices for FT from MAP. 
    tday.new.inds = MAP(remapped.new,MAPcols(2)); 
    tday.newtuning.inds = MAP(remapped.tuning,MAPcols(2));

%% Find peaks in tuning curve. 
    %Peaks in tuning curve (bins). 
    [~,tday.new.binpks] = cellfun(@max,CURVES{2}.tuning(tday.new.inds)); 
    [~,tday.newtuning.binpks] = cellfun(@max,CURVES{2}.tuning(tday.newtuning.inds)); 
    
    [~,yday.binpks] = cellfun(@max,CURVES{1}.tuning(TIMECELLS{1})); 
    
    %Peaks in tuning curve (time). 
    tday.new.tpks = t(tday.new.binpks);
    tday.newtuning.tpks = t(tday.newtuning.binpks); 
    
    yday.tpks = yt(yday.binpks); 
    
    %Scaling for the first plot. This is so scatter will plot circles with
    %diameter sized according to time peaks instead of area. 
    tday.new.tpks_scaled = 2*sqrt((tday.new.tpks+0.1)/pi)*30; 
    tday.newtuning.tpks_scaled = 2*sqrt((tday.newtuning.tpks+0.1)/pi)*30;
    
    yday.tpks_scaled = 2*sqrt((yday.tpks+0.1)/pi)*30; 

%% Plot some topography. 
    %Get the tday.centroids of the neurons. 
    tday.centroids = getNeuronCentroids(tday.Animal,tday.Date,tday.Session);
    yday.centroids = getNeuronCentroids(yday.Animal,yday.Date,yday.Session);
        
    load(fullfile(tday.Location,'ProcOut.mat'),'NeuronImage','Xdim','Ydim');
    
    %Plot neurons sized by when their time field occurs. Larger = later in
    %the delay. 
    figure;
    scatter(tday.centroids(tday.new.inds,1),...
        tday.centroids(tday.new.inds,2),...
        tday.new.tpks_scaled); hold on;
    scatter(tday.centroids(tday.newtuning.inds,1),...
        tday.centroids(tday.newtuning.inds,2),...
        tday.newtuning.tpks_scaled); hold on;
    scatter(yday.centroids(TIMECELLS{1},1),...
        yday.centroids(TIMECELLS{1},2),...
        yday.tpks_scaled); 
    xlim([0 Xdim]); ylim([0 Ydim]); axis off; 
    title('New Time Fields'); legend({'New Time Cells','New Tuning'}); 

    %Plot where the new time cells are. Red = new time cells. Blue = time
    %cells that changed their tuning. 
%     figure('position',[-1500 -40 1080 900]); 
%     PlotNeurons(NeuronImage,Xdim,Ydim,'g',0.5); hold on;
%     PlotNeurons(NeuronImage(tday.new.inds),Xdim,Ydim,'r',2); hold on;
%     PlotNeurons(NeuronImage(tday.newtuning.inds),Xdim,Ydim,'b',2);
%     title('New Time Tuning Topography'); 
    
%% 
    window = 0.1; 
    keepgoing = 1;
    i = 1;
    f = figure;
    while keepgoing
        yInfield = TIMECELLS{1}(yday.tpks < tday.new.tpks(i) + window & ...
            yday.tpks > tday.new.tpks(i) - window);
        
        scatter(tday.centroids(tday.new.inds(i),1),...
            tday.centroids(tday.new.inds(1),2),...
            tday.new.tpks_scaled(i)); hold on;
        scatter(yday.centroids(yInfield,1),...
            yday.centroids(yInfield,2),...
            yday.tpks_scaled(ismember(TIMECELLS{1},yInfield)));  
        xlim([0 Xdim]); ylim([0 Ydim]); axis off; 
        hold off; 
        
        [keepgoing,i] = scroll(i,length(tday.new.inds),f); 
        
    end
            
    i = 1;
    keepgoing = 1; 
    while keepgoing
        yInfield = TIMECELLS{1}(yday.tpks < tday.newtuning.tpks(i) + window & ...
            yday.tpks > tday.newtuning.tpks(i) - window);
        
        scatter(tday.centroids(tday.newtuning.inds(i),1),...
            tday.centroids(tday.newtuning.inds(1),2),...
            tday.newtuning.tpks_scaled(i)); hold on;
        scatter(yday.centroids(yInfield,1),...
            yday.centroids(yInfield,2),...
            yday.tpks_scaled(ismember(TIMECELLS{1},yInfield)));  
        xlim([0 Xdim]); ylim([0 Ydim]); axis off; 
        hold off; 
        
        [keepgoing,i] = scroll(i,length(tday.newtuning.inds),f); 
        
    end   
    keyboard; 
    
end