function [stability,mds] = stabilityMetric(mapMD,baseMD,compMDs,type)
%
%
%

%% Set up
    %Compile necessary data across sessions. 
    DATA = CompileMultiSessionData([baseMD,compMDs],{'t','timecells'});
    
    %Make sure all treadmill durations are equal across sessions. 
    Ts = cell2mat(DATA.t); 
    T = unique(Ts);
    assert(length(T)==1,'More than one unique run duration!');
    
    %Preallocate and useful variables. 
    mds = [baseMD,compMDs];
    dates = {mds.Date};
    
    %Sort chronologically.
    [~,order] = sort(datenum(dates,'mm_dd_yyyy')); 
    mds = mds(order); 
    
    nSessions = length(mds); 
    %Using the distance metric, find the distance between the peak of each
    %time cell to itself on other days. 
    if strcmp(type,'distance')
        [~,sortedPeaks] = multiPastalkovaPlot(mapMD,baseMD,compMDs,Ts,false); 

        nTCs = size(sortedPeaks,1);
        stability = nan(nTCs,nSessions); 
        stability(:,1) = zeros(nTCs,1); 
        %For each session, take absolute difference. 
        for s = 2:nSessions
            stability(:,s) = abs(sortedPeaks(:,s) - sortedPeaks(:,1));
        end
        
        %Reorder the columns so that they match the order in mds. 
        stability = stability(:,order); 
    else
        
    end
    
end