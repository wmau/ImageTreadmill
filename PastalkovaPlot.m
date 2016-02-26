function [sortedPastalkova,order] = PastalkovaPlot(animal,date,session,T,plotit)
%sortedPastalKova = PastalkovaPlot(animal,date,session,T)
%
%   Makes a plot that shows response curves that tile the delay. 

%% Make a plot that shows response curves that tile the delay. 
    ChangeDirectory(animal,date,session);

    try
        load('TimeCells.mat'); 
    catch
        [TimeCells,ratebylap,curves] = FindTimeCells(animal,date,session,T); 
    end 
    
    %Concatenate all time cell tuning curves.
    tilemat = cell2mat(curves.tuning(TimeCells));
    nBins = size(tilemat,2); 
    
    %Find the peak and normalize.
    [peaks,peakInds] = max(tilemat,[],2);
    normtilemat = tilemat./repmat(peaks,1,nBins);
    [~,order] = sort(peakInds);
    sortedPastalkova = normtilemat(order,:);
    
    %Plot. 
    if plotit
    figure;
    imagesc([0:T],[1:5:length(TimeCells)],sortedPastalkova);
        colormap gray;
        xlabel('Time [s]'); ylabel('Neurons'); 
    end

end