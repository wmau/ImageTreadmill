%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309); 
    
    cellType = 'time';
    infoType = 'ti';
    disregardStability = false;
    
    switch infoType
        case 'ti', infoTypeStr = 'Temporal';
        case 'si', infoTypeStr = 'Spatial';
    end
    
    switch cellType
        case 'time', c = [0 .5 .5]; otherC = [.58 .44 .86];
        case 'place',c = [.58 .44 .86]; otherC = [0 .5 .5];
        case 'dual',c = [0 0 0]; otherC = [.5 .5 .5];
    end
    
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    
    [alignedDiffs,stable,unstable,dayDelta,stableMeans,unstableMeans,alignedDays,alignedOtherI] = ...
        deal(cell(nAnimals,1)); 
    [I,days,stabilityLabel,otherI,otherIDays] = deal([]);
    figure; hold on;
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        %Align all the informations to peak.
        [alignedDiffs{a},stable{a},unstable{a},dayDelta{a},alignedDays{a},alignedOtherI{a}] = ...
            DiffInfoFromPeak(fulldataset(ssns),cellType,infoType,'plotit',false);
        
        I = [I; alignedDiffs{a}(:)];
        days = [days; alignedDays{a}(:)];
        
        otherI = [otherI; alignedOtherI{a}(:)];
        otherIDays = [otherIDays; alinedDays{a}(:)];
        
    end
    
    %Eliminate bad days.
    bad = isnan(days);
    days(bad) = [];
    I(bad) = [];
    
    %Do the same for the opposite dimension.
    bad = isnan(otherIDays); 
    otherIDays(bad) = [];
    otherI(bad) = [];
    
    %X-axis
    daysInd = days + abs(min(days)) + 1;
    otherIDaysInd = otherIDays + abs(min(otherIDays)) + 1; 
    dayDeltaUse = [min(days):max(days)];
    dayDeltaUseOther = [min(otherIDays):max(otherIDays)];
    
    %Get means and sem. 
    m = accumarray(daysInd,I,[],@nanmean);
    sem = accumarray(daysInd,I,[],@nanstd); 
    t = tabulate(daysInd); 
    sem = sem ./ sqrt(t(:,2)); 
    
    otherIMean = accumarray(daysInd,otherI,[],@nanmean);
    otherISEM = accumarray(daysInd,otherI,[],@nanstd); 
    t = tabulate(otherDaysInd); 
    otherISEM = otherISEM ./ sqrt(t(:,2)); 
    
    errorbar(dayDeltaUse+0.1,m,sem,'color',c,'linewidth',4);
    errorbar(dayDeltaUseOther-.1,otherIMean,otherISEM,'color',otherC,'linewidth',4);
    