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
    treadmillTime = 'early';
    
%% Step 1: Depict temporal information based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','TI',treadmillTime);
%     [S,US] = ParseInfoStabilityParams(fulldataset,'time','FR');
%     
%     figure;
%     scatter([s,us]',[S,US]',2,'filled');
%     [~,p] = corr([s,us]',[S,US]')
    
    fPos = [-1900 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});

%% Step 2: Depict temporal information based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','TI',treadmillTime);
    
    fPos = [-1600 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    
%% Step 3: Depict spatial information based on spatial stability. 
     [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','SI',treadmillTime);
%     [S,US] = ParseInfoStabilityParams(fulldataset,'place','FR');
%     
%     figure;
%     scatter([s,us]',[S,US]',2,'filled');
%     [~,p] = corr([s,us]',[S,US]')

    fPos = [-1300 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    
%% Step 4: Depict spatial information based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','SI',treadmillTime);

    fPos = [-1000 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    
%% Step 5: Depict firing rate based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','FR',treadmillTime);

    fPos = [-700 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',teal,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    
%% Step 6: Depict firing rate based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','FR',treadmillTime);

    fPos = [-400 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen''s d = ',num2str(d)]});
    
%%
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','pk',treadmillTime);

    fPos = [-400 460 300 450];
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Time Peak [s]','boxColor',purple,'position',...
        fPos,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    tp = ranksum(s,us); 
    d = cohensD(s,us);
    title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)], ['Cohen"s d = ',num2str(d)]});
    
%% 
function [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(data,stability,infoType,treadmillTime)
    nAnimals = length(unique({data.Animal}));
    colors = parula(nAnimals);
        
    [I,N] = StabilityStatsEarlyLate(data,stability,infoType,treadmillTime);
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