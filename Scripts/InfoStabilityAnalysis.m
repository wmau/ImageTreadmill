%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = [MD(292:303) MD(305:308)];   
    %fulldataset(9:13) = [];

    %Some initial variables. 
    teal = [0 .5 .5];
    purple = [.58 .44 .86];
    
%% Save information
    saveBool = false;
    folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures';
    TimeTI = fullfile(folder,'Stable Time TI');
    TimeSI = fullfile(folder,'Stable Time SI');
    TimeFR = fullfile(folder,'Stable Time FR');
    PlaceSI = fullfile(folder,'Stable Place SI');
    PlaceTI = fullfile(folder,'Stable Place TI');
    PlaceFR = fullfile(folder,'Stable Place FR'); 
    
    if saveBool
        c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');
    
        if ~strcmp(c,'y')
            saveBool = false;
        end
    end
    
%% Step 1: Depict temporal information based on temporal stability. 
    %Categorize.
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','TI');
%     [S,US] = ParseInfoStabilityParams(fulldataset,'time','FR');
%     
%     figure;
%     scatter([s,us]',[S,US]',2,'filled');
%     [~,p] = corr([s,us]',[S,US]')
    
    fPos = [-1900 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeTI,'-dpdf');
    end
    
%% Step 2: Depict spatial information based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','SI');

    fPos = [-1600 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeSI,'-dpdf');
    end    
    
%% Step 3: Depict firing rate based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','FR');

    fPos = [-1300 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(TimeFR,'-dpdf');
    end  
    
%% Step 4: Depict spatial information based on spatial stability. 
     [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','SI');
%     [S,US] = ParseInfoStabilityParams(fulldataset,'place','FR');
%     
%     figure;
%     scatter([s,us]',[S,US]',2,'filled');
%     [~,p] = corr([s,us]',[S,US]')

    fPos = [-1000 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceSI,'-dpdf');
    end  
  
%% Step 5: Depict temporal information based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','TI');
    
    fPos = [-700 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceTI,'-dpdf');
    end  
     
%% Step 6: Depict firing rate based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','FR');

    fPos = [-400 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    if saveBool
        print(PlaceFR,'-dpdf');
    end 
    
%%
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','fluor');

    fPos = [-400 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Fluor. Intensity [AU]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen"s d = ',num2str(d)]});

%% 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','fluor');

    fPos = [-400 460 300 450];
    scatterBox([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Norm. Fluor. Intensity [AU]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen"s d = ',num2str(d)]});

%% Nested function
function [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(data,stability,infoType)
    nAnimals = length(unique({data.Animal}));
    colors = parula(nAnimals);
        
    [I,N] = PartitionStats(data,stability,infoType);
    s = cell2mat(I.stable')';          %ALL stable stats.
    us = cell2mat(I.unstable')';       %ALL unstable stats.
    sN = length(s);                    %Number of stable neurons.
    usN = length(us);                  %Number of unstable neurons. 

    %Group identity vector. 
    grps = [zeros(1,sN), ones(1,usN)];
    sColors = []; 
    usColors = [];

    for a=1:nAnimals     
        sColors = [sColors; repmat(colors(a,:),N.stable(a),1)];
        usColors = [usColors; repmat(colors(a,:),N.unstable(a),1)];
    end
    
end