function toptempdists(animal,date,session,T)
%function toptempdists(animal,date,session,T)
%
%   Rudimentary version of a function to see whether there is a correlation
%   between how far apart time cell peaks are and their physical distances
%   in the hippocampus. 

%%
    ChangeDirectory(animal,date,session);

    try
        load('TimeCells.mat'); 
    catch
        [TimeCells,ratebylap,curves,delays,x,y,time_interp] = FindTimeCells(animal,date,session,T); 
    end
    
    centroids = getNeuronCentroids(animal,date,session);
    
    tempResponses = cell2mat(curves.tuning(TimeCells)); 
    [~,peakT] = max(tempResponses,[],2); 
    
    tempdist = pdist([peakT, zeros(length(peakT),1)])/20;
    topdist = pdist(centroids(TimeCells,:));
    
    scatter(tempdist,topdist,'.');
    [r,p] = corr(tempdist',topdist');
    
    r
    p
end