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
    fPos = [-1300 -55 620 920];
    figure('Position',fPos);
    teal = [0 .5 .5];
    purple = [.58 .44 .86];
    
%% Step 1: Depict temporal information based on temporal stability. 
    %Categorize.
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','TI');
    
    subplot(3,2,1);
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',teal,'position',...
        false,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    %tp = ranksum(s,us); 
    title(['KS p = ',num2str(kp)]);
    
%% Step 2: Depict temporal information based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','TI');
    
    subplot(3,2,2);
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Temporal Information [bits]','boxColor',purple,'position',...
        false,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    %tp = ranksum(s,us); 
    title(['KS p = ',num2str(kp)]);
    
%% Step 3: Depict spatial information based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','SI');

    subplot(3,2,3);
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',purple,'position',...
        false,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    %tp = ranksum(s,us); 
    title(['KS p = ',num2str(kp)]);
    
%% Step 4: Depict spatial information based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','SI');

    subplot(3,2,4);
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Spatial Information [bits]','boxColor',teal,'position',...
        false,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    %tp = ranksum(s,us); 
    title(['KS p = ',num2str(kp)]);
    
%% Step 5: Depict firing rate based on temporal stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'time','FR');

    subplot(3,2,5);
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',teal,'position',...
        false,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    %tp = ranksum(s,us); 
    title(['KS p = ',num2str(kp)]);
    
%% Step 6: Depict firing rate based on spatial stability. 
    [s,us,grps,sColors,usColors] = ParseInfoStabilityParams(fulldataset,'place','FR');

    subplot(3,2,6);
    boxScatterplot([s,us],grps,'xLabels',{'Stable','Unstable'},...
        'yLabel','Transient Frequency','boxColor',purple,'position',...
        false,'circleColors',[sColors;usColors]);
    [~,kp] = kstest2(s,us);
    %tp = ranksum(s,us); 
    title(['KS p = ',num2str(kp)]);
    
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