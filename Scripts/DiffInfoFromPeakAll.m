%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = [MD(292:299) MD(300:303) MD(305:308)];   
    
    cellType = 'time';
    infoType = 'ti';
    disregardStability = true;
    
    switch infoType
        case 'ti', infoTypeStr = 'Temporal';
        case 'si', infoTypeStr = 'Spatial';
    end
    
    switch cellType
        case 'time', c = [0 .5 .5]; otherC = [.58 .44 .86];
        case 'place',c = [.58 .44 .86]; otherC = [0 .5 .5];
        case {'all','dual'},c = [0 0 0]; otherC = [.5 .5 .5];
    end
    
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    
    [alignedDiffs,stable,unstable,dayDelta,stableMeans,unstableMeans,alignedDays,alignedOtherI] = ...
        deal(cell(nAnimals,1)); 
    [I,days,stabilityLabel,otherI,otherIDays] = deal([]);
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        %Align all the informations to peak.
        [alignedDiffs{a},stable{a},unstable{a},dayDelta{a},alignedDays{a},alignedOtherI{a}] = ...
            DiffInfoFromPeak(fulldataset(ssns),cellType,infoType,'plotit',false);
        
        %Grab stable information and day vectors. 
        tempI = alignedDiffs{a}(stable{a},:); 
        tempI = tempI(:);
        tempDays = alignedDays{a}(stable{a},:);
        tempDays = tempDays(:); 
             
        I = [I; tempI];
        days = [days; tempDays];
        stabilityLabel = [stabilityLabel; ones(size(tempI))];
        
        %Grab unstable. 
        tempI = alignedDiffs{a}(unstable{a},:);
        tempI = tempI(:);
        tempDays = alignedDays{a}(unstable{a},:);
        tempDays = tempDays(:); 
        
        I = [I; tempI];
        days = [days; tempDays];
        stabilityLabel = [stabilityLabel; zeros(size(tempI))];
        
        %Get other dimension of information.
        tempOtherI = alignedOtherI{a};
        tempOtherI = tempOtherI(:);
        tempDays = alignedDays{a};
        tempDays = tempDays(:);
        otherIDays = [otherIDays; tempDays]; 
        otherI = [otherI; tempOtherI];
    end
    
    %Eliminate bad days.
    bad = isnan(days);
    days(bad) = [];
    I(bad) = [];
    stabilityLabel(bad) = [];
    
    %Do the same for the opposite dimension.
    bad = isnan(otherIDays); 
    otherIDays(bad) = [];
    otherI(bad) = [];
    
    daysInd = days + abs(min(days)) + 1;
    otherIDaysInd = otherIDays + abs(min(otherIDays)) + 1; 
    stable = stabilityLabel == 1;
    unstable = stabilityLabel == 0;
    dayDeltaUse = [min(days):max(days)];
    dayDeltaUseOther = [min(otherIDays):max(otherIDays)];
    
    %Get means and SEM.
    stableMean = accumarray(daysInd(stable),I(stable),[],@nanmean);
    stableSEM = accumarray(daysInd(stable),I(stable),[],@nanstd);
    t = tabulate(daysInd(stable));
    stableSEM = stableSEM ./sqrt(t(:,2)); 
  
    unstableMean = accumarray(daysInd(unstable),I(unstable),[],@nanmean);
    unstableSEM = accumarray(daysInd(unstable),I(unstable),[],@nanstd); 
    t = tabulate(daysInd(unstable));
    unstableSEM = unstableSEM ./sqrt(t(:,2)); 
    
    otherIMean = accumarray(otherIDaysInd,otherI,[],@nanmean);
    otherISEM = accumarray(otherIDaysInd,otherI,[],@nanstd); 
    t = tabulate(otherIDaysInd); 
    otherISEM = otherISEM ./ sqrt(t(:,2)); 

    %Plot means across animals.
    figure; hold on;
    errorbar(dayDeltaUseOther-.1,otherIMean,otherISEM,'color',otherC,'linewidth',4,...
        'capsize',0);
    set(gca,'tickdir','out','fontsize',12,'linewidth',4); 
    ylabel('Fraction of Peak Info.','fontsize',15);
    xlabel('Days from Peak','fontsize',15);
    
    if disregardStability
        allDays = cell2mat(cellfun(@(x) x(:),alignedDays,'unif',0));
        allIs = cell2mat(cellfun(@(x) x(:),alignedDiffs,'unif',0));
        
        dayX = [min(allDays):max(allDays)];
        allDaysInd = allDays + abs(min(allDays)) + 1;
        bad = isnan(allDaysInd); 
        allDaysInd(bad) = [];
        allIs(bad) = [];
        
        m = accumarray(allDaysInd,allIs,[],@nanmean);
        sem = accumarray(allDaysInd,allIs,[],@nanstd); 
        t = tabulate(allDays);
        sem = sem./sqrt(t(:,2)); 
              
        errorbar(dayX+.1,m,sem,'color',c,'linewidth',4,'capsize',0);
        
        withinDim = [ones(size(I)); zeros(size(otherI))];
        DAYS = [days; otherIDays];
        Is = [I; otherI];
        bad = DAYS==-4 & withinDim==0;
        DAYS(bad) = []; withinDim(bad) = []; 
        Is(bad) = [];
        
        grps = {DAYS, withinDim};
        
        [~,tbl,stats] = anovan(Is,grps,'model','full','varnames',{'Days','Within Dimension'});
        figure;
        multcompare(stats,'dimension',[1,2],'ctype','scheffe');
    else
        errorbar(dayDeltaUse+.1,stableMean,stableSEM,'color',c,'linewidth',4,'capsize',0);
        errorbar(dayDeltaUse,unstableMean,unstableSEM,'color','r','linewidth',4,'capsize',0);
    
        %ANOVA.
        grps = {days,stabilityLabel};
        [~,tbl,stats] = anovan(I,grps,'model','full','varnames',{'Day','Stability'});
        comps = multcompare(stats,'dimension',[1,2],'display','off','ctype','scheffe');
        disp(['Stable vs. unstable on Day -3 P = ',num2str(comps(8,6))]);
        disp(['Stable vs. unstable on Day -2 P = ',num2str(comps(23,6))]);
        disp(['Stable vs. unstable on Day -1 P = ',num2str(comps(37,6))]);
        disp(['Stable vs. unstable on Day 0 P = ',num2str(comps(50,6))]);
        disp(['Stable vs. unstable on Day +1 P = ',num2str(comps(62,6))]);
        disp(['Stable vs. unstable on Day +2 P = ',num2str(comps(73,6))]);
        disp(['Stable vs. unstable on Day +3 P = ',num2str(comps(83,6))]);
        disp(['Stable vs. unstable on Day +4 P = ',num2str(comps(92,6))]);
    
    end
    
    
