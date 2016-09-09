function PlotStability(mds)
%
%
%

%% 
    %What animals are in md?
    animals = unique({mds.Animal}); 
    nAnimals = length(animals); 
    
    DATA = CompileMultiSessionData(mds,{'t'});
    
    %Get mapMD for each animal.
    for a=1:nAnimals
        %Get all the MDs for that animal and their corresponding mapMDs. 
        mdIdx = ismember({mds.Animal},animals{a});
        theseMDs = mds(mdIdx);
        mapMDs(a) = getMapMD(theseMDs);          
    end
    
    %Preallocate things. 
    stability =     cell(1,nAnimals); 
    meanStability = cell(1,nAnimals); 
    dayrank =       cell(1,nAnimals); 
    STABILITY = [];
    DAYRANK = []; 
    
    figure; 
    for a=1:nAnimals
        %Get all the MDs for that animal.
        mdIdx = ismember({mds.Animal},animals{a});
        theseMDs = mds(mdIdx);
        
        %Make sure the dates are sorted chronologically. 
        dates = {theseMDs.Date};
        theseMDs = sortbyDate(dates,theseMDs); 
        
        %Get the middle MDs and the rest. 
        middle = ceil(length(theseMDs)/2);
        baseMD = theseMDs(middle);
        compMDs = theseMDs([1:middle-1,middle+1:end]);
        
        %Get neurons across days and sort according to the base session.
        Ts = [DATA.t{mdIdx}];
        [~,sortedPeaks] = multiPastalkovaPlot(mapMDs(a),baseMD,compMDs,Ts,false);
        
        %Compute stability. 
        stability{a} = sDistanceMetric(sortedPeaks);
        stability{a} = sortbyDate({baseMD.Date compMDs.Date},stability{a});
        STABILITY = [STABILITY; stability{a}(:)];   %Linearize for accumarray later. 
        
        %Take the mean. 
        meanStability{a} = nanmean(stability{a}); 
        
        %Linearize for accumarray. 
        dayrank{a} = zeros(size(stability{a}));     
        for s=1:length(theseMDs)
            dayrank{a}(:,s) = daysact(baseMD.Date,theseMDs(s).Date); 
        end
        DAYRANK = [DAYRANK; dayrank{a}(:)];
        
        %Plot animal points. 
        scatter(unique(dayrank{a}),meanStability{a},'x'); hold on; 
    end
 
    %Referencing the values in STABILITY, the days in chronological order.
    %This line makes these values non-negative. 
    DAYRANK = DAYRANK + abs(min(DAYRANK))+1;
    
    %Take the mean for each session. 
    MEANSTABILITY = accumarray(DAYRANK,STABILITY,[],@nanmean);
    
    %X-axis from -X to +Y days from reference. 
    XAXISDATE = min(cellfun(@min,cellfun(@min,dayrank,'unif',false)))...
        :max(cellfun(@max,cellfun(@max,dayrank,'unif',false)));
    
    %Plot. 
    plot(XAXISDATE,MEANSTABILITY); 
    keyboard;
    
end