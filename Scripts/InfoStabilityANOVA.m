%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309);      

    %Some initial variables. 
    teal = [0 .5 .5];
    purple = [.58 .44 .86];
    
%% Save information
    saveBool = true;
    folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\';
    Time = fullfile(folder,'Stable Time');
    TimeFR = fullfile(folder,'Stable Time FR');
    Place = fullfile(folder,'Stable Place');
    PlaceFR = fullfile(folder,'Stable Place FR'); 
    
    if saveBool
        c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
        if ~strcmp(c,'y')
            saveBool = false;
        end
    end
    
    %Main analysis.
    figure('Position',[130 240 250 440]);  hold on;
    [sSame,usSame] = ParseInfoStability(fulldataset,'time','ti',teal);
    [sDiff,usDiff] = ParseInfoStability(fulldataset,'time','si',purple);
    MakeErrBar(sDiff,usDiff,purple);
    MakeErrBar(sSame,usSame,teal);
    ylabel('z-scored Mutual Info.');
    xlabel('Temporal Stability','Color',teal);
    if saveBool
        print(Time,'-dpdf');
    end
    
    %Run ANOVA.
    disp('Running Two-Way ANOVA of MI based on temporal stability.');
    ANOVAit(sSame,usSame,sDiff,usDiff);
    
    %Check calcium event frequency.
    figure('Position',[400 240 250 440]); hold on;
    [S,US] = ParseInfoStability(fulldataset,'time','fr',[.7 .7 .7]);
    MakeErrBar(S,US,[.7 .7 .7]);
    title(['P = ',num2str(ranksum(S,US))]);
    ylabel('Norm. Ca Event Freq.');
    xlabel('Temporal Stability','Color',teal);
    if saveBool
        print(TimeFR,'-dpdf');
    end
    
    
%% For space.
    disp(' ');
    
    %Main analysis.
    figure('Position',[670 240 250 440]); hold on;
    [sSame,usSame] = ParseInfoStability(fulldataset,'place','si',purple);
    [sDiff,usDiff] = ParseInfoStability(fulldataset,'place','ti',teal);   
    MakeErrBar(sDiff,usDiff,teal);
    MakeErrBar(sSame,usSame,purple);
    ylabel('z-scored Mutual Info.');
    xlabel('Spatial Stability','Color',purple);
    if saveBool
        print(Place,'-dpdf');
    end
    
    %Run ANOVA.
    disp('Running Two-Way ANOVA of MI based on spatial stability.');
    ANOVAit(sSame,usSame,sDiff,usDiff);
    
    %Check calcium event frequency.
    figure('Position',[940 240 250 440]); hold on;
    [S,US] = ParseInfoStability(fulldataset,'place','fr',[.7 .7 .7]);
    
    MakeErrBar(S,US,[.7 .7 .7]);
    title(['P = ',num2str(ranksum(S,US))]);
    ylabel('Norm. Ca Event Freq.');
    xlabel('Spatial Stability','Color',purple);
    if saveBool
        print(PlaceFR,'-dpdf');
    end

%%
function [sM,usM] = ParseInfoStability(fulldataset,stabilityType,infoType,c)
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);

    [stats,~,stable,unstable] = PartitionStats(fulldataset,stabilityType,infoType);
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        [i,j] = deal(1);
        for s=1:length(ssns)-1
            sExtent = length(stable{a}{s}); 
            usExtent = length(unstable{a}{s});
            
            if sExtent > 1, sM{a}(s) = median(stats.stable{a}(i:i+sExtent-1));
            else, sM{a}(s) = nan; end
                
            if usExtent > 1, usM{a}(s) = median(stats.unstable{a}(j:j+usExtent-1)); 
            else, usM{a}(s) = nan; end
            
            i = i+sExtent;
            j = j+usExtent;
        end
    end
    
    sM = cell2mat(sM)';
    usM = cell2mat(usM)';
    p = signrank(sM,usM);  
    
    plot([1,2],[sM,usM],'color',[c .5]);
    set(gca,'xticklabel',{'Stable','Unstable'});
    set(gca,'xtick',[1:2],'tickdir','out');
    %title(['P = ',num2str(p)]);
    xlim([0.5,2.5]);
end

%%
function [sM,usM,sSEM,usSEM] = MakeErrBar(stable,unstable,c)
    nSessions = length(stable); 
    sM = nanmean(stable);          sSEM  = nanstd(stable)./sqrt(nSessions);
    usM = nanmean(unstable);       usSEM = nanstd(unstable)./sqrt(nSessions); 
    
    errorbar([1,2],[sM,usM],[sSEM,usSEM],'color',c,'linewidth',5);
end

%%
function ANOVAit(sSame,usSame,sDiff,usDiff)
    nSessions = length(sSame);
    
    dimension = [ones(nSessions*2,1); zeros(nSessions*2,1)];
    stability = repmat([ones(nSessions,1); zeros(nSessions,1)],[2,1]);
    grps = {dimension,stability};
    X = [sSame; usSame; sDiff; usDiff];
    [~,tbl,stats] = anovan(X,grps,'model','full','varnames',...
        {'Dimension','Stability'},'display','off');
    comps = multcompare(stats,'dimension',[1,2],'display','off','ctype','hsd');
    
    disp(['Main effect of dimension F = ',num2str(tbl{2,6}),', P = ',num2str(tbl{2,7})]);
    disp(['Main effect of stability F = ',num2str(tbl{3,6}),', P = ',num2str(tbl{3,7})]);
    disp(['Interaction F = ',num2str(tbl{4,6}),', P = ',num2str(tbl{4,7})]);
    disp(['Difference between stable and unstable, same dimension P = ',num2str(comps(5,6))]);
    disp(['Difference between stable and unstable, opposite dimension P = ',num2str(comps(2,6))]);
    
end