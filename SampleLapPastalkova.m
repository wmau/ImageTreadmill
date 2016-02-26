function SampleLapPastalkova(animal,date,session,T) 
%SampleLapPastalkova(animal,date,session,T) 
%
%   Per Nick's suggestion, how do we know the tuning curves are reliable?
%   For one, I have a consistency filter (at the time of this writing, must
%   fire on 1/5 of all laps). This is another demonstration of time cell
%   reliabilty. Randomly sample laps and calculate new tuning curves (trial
%   means) for the partial datasets. Then make Pastalkova plot of the
%   complete dataset and compare to the randomly sampled ones, plotted in
%   the same order. If time cells are reliable, diagonal should be
%   conserved. 
%
%   INPUTS: 
%       animal,date,session,T: Which session to plot. 
%

%% Get the Pastalkova plot for the entirety of neurons. 
    [sortedPastalkova,order] = PastalkovaPlot(animal,date,session,T,0); 

%% Randomly sample laps. 
    %Load time cell indices, trial by trial rate maps, and data about delay
    %duration and complete runs. 
    load(fullfile(pwd,'TimeCells.mat'),'TimeCells','ratebylap','TodayTreadmillLog'); 
    
    %Find all lap numbers that are valid. 
    lapNumComplete = find(TodayTreadmillLog.complete & TodayTreadmillLog.delaysetting==T); 
    nComplete = length(lapNumComplete); 
    
    %Preallocate. 
    even = logical(zeros(nComplete,1)); 
    
    %Randomly sample good laps. 
    even(randsample(nComplete,round(nComplete/2))) = true; 
    odd = ~even;
    
    %Get lap numbers. 
    evenLaps = lapNumComplete(even); 
    oddLaps = lapNumComplete(odd); 

%% Get tuning curves from partial data set and plot. 
    %Take mean of select laps. 
    eTuning = squeeze(mean(ratebylap(evenLaps,:,:)))';
    oTuning = squeeze(mean(ratebylap(oddLaps,:,:)))'; 
    
    %Number of time bins. 
    nBins = size(eTuning,2); 
    
    %Normalization. 
    ePeaks = max(eTuning,[],2); 
    oPeaks = max(oTuning,[],2); 
    
    eTuning = eTuning./repmat(ePeaks,1,nBins); 
    oTuning = oTuning./repmat(oPeaks,1,nBins); 

    %Plot. 
    figure;
    subplot(1,3,1);
        imagesc([0:T],[1:5:length(TimeCells)],sortedPastalkova);
            colormap gray; ylabel('Neurons'); 
    %"Even/odd" trials. 
    subplot(1,3,2); 
        imagesc([0:T],[1:5:length(TimeCells)],eTuning(TimeCells(order),:));
            colormap gray; title('Even Trials'); xlabel('Time [s]');
    subplot(1,3,3); 
        imagesc([0:T],[1:5:length(TimeCells)],oTuning(TimeCells(order),:));
            colormap gray; title('Odd Trials'); 

end