%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309); 
    
    cellType = 'place';
    infoType = 'si';
    
    switch cellType
        case 'time', c = [0 .5 .5];
        case 'place',c = [.58 .44 .86];
    end
    
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    
    [alignedDiffs,stable,unstable,dayDelta,stableMeans,unstableMeans,alignedDays] = ...
        deal(cell(nAnimals,1)); 
    [I,days,stabilityLabel] = deal([]);
    figure; hold on;
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        %Align all the informations to peak.
        [alignedDiffs{a},stable{a},unstable{a},dayDelta{a},alignedDays{a}] = ...
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

    end
    
    %Eliminate bad days.
    bad = isnan(days);
    days(bad) = [];
    daysInd = days + abs(min(days)) + 1;
    I(bad) = [];
    stabilityLabel(bad) = [];
    stable = stabilityLabel == 1;
    unstable = stabilityLabel == 0;
    dayDeltaUse = [min(days):max(days)];
    
    %Get means and SEM.
    stableMean = accumarray(daysInd(stable),I(stable),[],@nanmean);
    stableSEM = accumarray(daysInd(stable),I(stable),[],@nanstd);
    t = tabulate(daysInd(stable));
    stableSEM = stableSEM ./sqrt(t(:,2)); 
  
    unstableMean = accumarray(daysInd(unstable),I(unstable),[],@nanmean);
    unstableSEM = accumarray(daysInd(unstable),I(unstable),[],@nanstd); 
    t = tabulate(daysInd(unstable));
    unstableSEM = unstableSEM ./sqrt(t(:,2)); 

    %Plot means across animals.
    errorbar(dayDeltaUse,stableMean,stableSEM,'color',c,'linewidth',4);
    errorbar(dayDeltaUse,unstableMean,unstableSEM,'color','r','linewidth',4);
    set(gca,'tickdir','out','fontsize',12,'linewidth',4); 
    xlabel('Days from Peak','fontsize',15);
    ylabel('Fraction of Peak Info.','fontsize',15);
    
    %ANOVA.
    grps = {days,stabilityLabel};
    [~,tbl,stats] = anovan(I,grps,'model','full','varnames',{'Day','Stability'});
    comps = multcompare(stats,'dimension',[1,2],'display','off');
    disp(['Stable vs. unstable on Day -3 P = ',num2str(comps(8,6))]);
    disp(['Stable vs. unstable on Day -2 P = ',num2str(comps(23,6))]);
    disp(['Stable vs. unstable on Day -1 P = ',num2str(comps(37,6))]);
    disp(['Stable vs. unstable on Day -0 P = ',num2str(comps(50,6))]);
    disp(['Stable vs. unstable on Day +1 P = ',num2str(comps(62,6))]);
    disp(['Stable vs. unstable on Day +2 P = ',num2str(comps(73,6))]);
    disp(['Stable vs. unstable on Day +3 P = ',num2str(comps(83,6))]);
    disp(['Stable vs. unstable on Day +4 P = ',num2str(comps(92,6))]);
