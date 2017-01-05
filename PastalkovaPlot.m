function [sortedPastalkova,order] = PastalkovaPlot(MD,plotit)
%sortedPastalKova = PastalkovaPlot(animal,date,session,T)
%
%   Makes a plot that shows response curves that tile the delay. 

%% Make a plot that shows response curves that tile the delay. 
    cd(MD.Location); 
    
    try
        load('TimeCells.mat','TimeCells','curves','T'); 
    catch
        [TimeCells,ratebylap,curves] = FindTimeCells(MD,10); 
        tempInfo(MD);        
    end 
    load('TemporalInfo.mat','sig');
    TimeCells = intersect(find(sig),TimeCells);
    
    %Concatenate all time cell tuning curves.
    tilemat = cell2mat(curves.tuning(TimeCells));
    nBins = size(tilemat,2); 
    
    %Find the peak and normalize.
    [peaks,peakInds] = max(tilemat,[],2);
    normtilemat = tilemat./repmat(peaks,1,nBins);
    
    %Sort. 
    [~,order] = sort(peakInds);
    sortedPastalkova = normtilemat(order,:);
    
    %Plot. 
    if plotit
    figure;
    imagesc([0:T],[1:length(TimeCells)],sortedPastalkova);
        colormap gray;
        xlabel('Time [s]'); ylabel('Neurons'); 
    end

end