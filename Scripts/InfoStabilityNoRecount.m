%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309);      
    %fulldataset(9:13) = [];

    %Some initial variables. 
    teal = [0 .5 .5];
    purple = [.58 .44 .86];
    
    saveBool = true;
    folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals\Alternative Info Stability Analyses';
    TimeTI = fullfile(folder,'Stable Time TI no recount');
    TimeSI = fullfile(folder,'Stable Time SI no recount');
    TimeFR = fullfile(folder,'Stable Time FR no recount');
    PlaceSI = fullfile(folder,'Stable Place SI no recount');
    PlaceTI = fullfile(folder,'Stable Place TI no recount');
    PlaceFR = fullfile(folder,'Stable Place FR no recount'); 
    
    if saveBool
        c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
        if ~strcmp(c,'y')
            saveBool = false;
        end
    end
    
%% 
    [STATS,stabilityStatus] = StabilityStats(fulldataset,'time','ti');
    [s,us,grps,sColors,usColors] = UnpackStats(STATS);
    
    fPos = [-1900 460 300 450];
    boxScatterplot([s;us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Temporal Information [z-bits]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeTI,'-dpdf');
    end
    
%% 
    [STATS,stabilityStatus] = StabilityStats(fulldataset,'time','si');
    [s,us,grps,sColors,usColors] = UnpackStats(STATS);
    
    fPos = [-1600 460 300 450];
    boxScatterplot([s;us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Spatial Information [z-bits]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeSI,'-dpdf');
    end
%% 
    [STATS,stabilityStatus] = StabilityStats(fulldataset,'time','fr');
    [s,us,grps,sColors,usColors] = UnpackStats(STATS);
    
    fPos = [-1300 460 300 450];
    boxScatterplot([s;us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Ca Event Freq.','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeFR,'-dpdf');
    end
    
%% 
    [STATS,stabilityStatus] = StabilityStats(fulldataset,'place','si');
    [s,us,grps,sColors,usColors] = UnpackStats(STATS);
    
    fPos = [-1000 460 300 450];
    boxScatterplot([s;us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Spatial Information [z-bits]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceSI,'-dpdf');
    end
    
%% 
    [STATS,stabilityStatus] = StabilityStats(fulldataset,'place','ti');
    [s,us,grps,sColors,usColors] = UnpackStats(STATS);
    
    fPos = [-700 460 300 450];
    boxScatterplot([s;us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Temporal Information [z-bits]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceTI,'-dpdf');
    end
    
%% Step 6: Depict firing rate based on spatial stability. 
    [STATS,stabilityStatus] = StabilityStats(fulldataset,'place','fr');
    [s,us,grps,sColors,usColors] = UnpackStats(STATS);

    fPos = [-400 460 300 450];
    boxScatterplot([s;us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Ca Event Freq.','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceFR,'-dpdf');
    end
%%
function [s,us,grps,sColors,usColors] = UnpackStats(STATS)
    nAnimals = length(STATS.stable);
    colors = parula(nAnimals);
    
    s = []; us = [];
    for a=1:nAnimals
        s = [s; cell2mat(STATS.stable{a}')];
        us = [us; cell2mat(STATS.unstable{a}')];
    end
    
    sN = length(s);
    usN = length(us);
    
    %Group identity vector. 
    grps = [zeros(1,sN), ones(1,usN)];
    sColors = []; 
    usColors = [];
    
    for a=1:nAnimals     
        nStable = sum(cellfun('length',STATS.stable{a}));
        nUnstable = sum(cellfun('length',STATS.unstable{a}));
        
        sColors = [sColors; repmat(colors(a,:),nStable,1)];
        usColors = [usColors; repmat(colors(a,:),nUnstable,1)];
    end
end
