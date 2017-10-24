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
    
%% Save information
    saveBool = true;
    folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures';
    TimeTI = fullfile(folder,'Stable Time TI Downsample');
    TimeSI = fullfile(folder,'Stable Time SI Downsample');
    TimeFR = fullfile(folder,'Stable Time FR Downsample');
    PlaceSI = fullfile(folder,'Stable Place SI Downsample');
    PlaceTI = fullfile(folder,'Stable Place TI Downsample');
    PlaceFR = fullfile(folder,'Stable Place FR Downsample'); 
    
    if saveBool
        c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
        if ~strcmp(c,'y')
            saveBool = false;
        end
    end
    
%% Step 1: Depict temporal information based on temporal stability. 
    %Categorize.
    nStable = 165;
    nUnstable = 81;
    
    [s,us,grps] = ParseInfoStabilityParams(fulldataset,'time','TI',nStable,nUnstable);
%     [S,US] = ParseInfoStabilityParams(fulldataset,'time','FR');
%     
%     figure;
%     scatter([s,us]',[S,US]',2,'filled');
%     [~,p] = corr([s,us]',[S,US]')
    
    fPos = [-1900 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',teal,'position',...
        fPos);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeTI,'-dpdf');
    end
    
%% Step 2: Depict spatial information based on temporal stability. 
    [s,us,grps] = ParseInfoStabilityParams(fulldataset,'time','SI',nStable,nUnstable);

    fPos = [-1600 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',teal,'position',...
        fPos);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeSI,'-dpdf');
    end    
    
%% Step 3: Depict firing rate based on temporal stability. 
    [s,us,grps] = ParseInfoStabilityParams(fulldataset,'time','FR',nStable,nUnstable);
    
    fPos = [-1300 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',teal,'position',...
        fPos);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeFR,'-dpdf');
    end  
    
%% Step 4: Depict spatial information based on spatial stability. 
     [s,us,grps] = ParseInfoStabilityParams(fulldataset,'place','SI',nStable,nUnstable);
%     [S,US] = ParseInfoStabilityParams(fulldataset,'place','FR');
%     
%     figure;
%     scatter([s,us]',[S,US]',2,'filled');
%     [~,p] = corr([s,us]',[S,US]')

    fPos = [-1000 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',purple,'position',...
        fPos);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceSI,'-dpdf');
    end  
  
%% Step 5: Depict temporal information based on spatial stability. 
    [s,us,grps] = ParseInfoStabilityParams(fulldataset,'place','TI',nStable,nUnstable);
    
    fPos = [-700 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',purple,'position',...
        fPos);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceTI,'-dpdf');
    end  
     
%% Step 6: Depict firing rate based on spatial stability. 
    [s,us,grps] = ParseInfoStabilityParams(fulldataset,'place','FR',nStable,nUnstable);

    fPos = [-400 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',purple,'position',...
        fPos);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceFR,'-dpdf');
    end 
    
%%
%     [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','pk');
% 
%     fPos = [-400 460 300 450];
%     scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
%         'yLabel','Time Peak [s]','boxColor',purple,'position',...
%         fPos,'circleColors',[sColors;usColors]);
%     [~,kp] = kstest2(s,us);
%     tp = ranksum(s,us); 
%     d = cohensD(s,us);
%     title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen"s d = ',num2str(d)]});
    
%% Nested function
function [s,us,grps] = ParseInfoStabilityParams(data,stability,infoType,sN,usN)
    nAnimals = length(unique({data.Animal}));
    colors = parula(nAnimals);
        
    I = PartitionStats(data,stability,infoType);
    s = cell2mat(I.stable')';          %ALL stable stats.
    us = cell2mat(I.unstable')';       %ALL unstable stats.
    
    s = randsample(s,sN);
    us = randsample(us,usN);
    
    %Group identity vector. 
    grps = [zeros(1,sN), ones(1,usN)];

end