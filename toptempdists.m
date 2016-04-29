function toptempdists(MD,T)
%function toptempdists(animal,date,session,T)
%
%   Rudimentary version of a function to see whether there is a correlation
%   between how far apart time cell peaks are and their physical distances
%   in the hippocampus. 

%%
    animal = MD.Animal;
    date = MD.Date;
    session = MD.Session;
    ChangeDirectory(animal,date,session);

    try
        load('TimeCells.mat'); 
    catch
        [TimeCells,ratebylap,curves,delays,x,y,time_interp] = FindTimeCells(MD,T); 
    end
    
    centroids = getNeuronCentroids(MD);
    
    tempResponses = cell2mat(curves.tuning(TimeCells)); 
    [~,peakT] = max(tempResponses,[],2); 
    
    tempdist = pdist([peakT, zeros(length(peakT),1)]);
    topdist = pdist(centroids(TimeCells,:));
    
    scatter(tempdist,topdist,'.');
    [r,p] = corr(tempdist',topdist');
    
    r
    p
    keyboard;
end