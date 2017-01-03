%Full data set.
fulldataset = [MD(292:309)];
nAnimals = length(unique({fulldataset.Animal})); 
colors = parula(nAnimals);

%Partition temporal information scores into stable vs unstable based on
%time field or place field correlation.
[stblTimeI,stblTimeN] = PartitionStats(fulldataset,'time','TI'); 
[stblPlaceI,stblPlaceN] = PartitionStats(fulldataset,'place','TI'); 

%Stable time fields.
sTI_t = cell2mat(stblTimeI.stable')';
usTI_t = cell2mat(stblTimeI.unstable')';
grps_t = [zeros(1,length(sTI_t)), ones(1,length(usTI_t))];
animalColors_t_stable = nan(length(sTI_t),3);
animalColors_t_unstable = nan(length(usTI_t),3); 
s = 1;
u = 1;
for a = 1:nAnimals
    animalColors_t_stable(s:s+stblTimeN.stable(a)-1,:) = repmat(colors(a,:),...
        stblTimeN.stable(a),1);
    animalColors_t_unstable(u:u+stblTimeN.unstable(a)-1,:) = repmat(colors(a,:),...
        stblTimeN.unstable(a),1); 
    
    s = s+stblTimeN.stable(a);
    u = u+stblTimeN.unstable(a);
end
    
%Stable place fields.
sTI_p = cell2mat(stblPlaceI.stable')';
usTI_p = cell2mat(stblPlaceI.unstable')';
grps_p = [zeros(1,length(sTI_p)), ones(1,length(usTI_p))];
animalColors_p_stable = nan(length(sTI_p),3);
animalColors_p_unstable = nan(length(usTI_p),3); 
s = 1; 
u = 1;
for a = 1:nAnimals
    animalColors_p_stable(s:s+stblPlaceN.stable(a)-1,:) = repmat(colors(a,:),...
        stblPlaceN.stable(a),1);
    animalColors_p_unstable(u:u+stblPlaceN.unstable(a)-1,:) = repmat(colors(a,:),...
        stblPlaceN.unstable(a),1); 
    
    s = s+stblPlaceN.stable(a);
    u = u+stblPlaceN.unstable(a);
end

%%
[stblPlaceSpatialInfo,stblPlaceSpatialN] = PartitionStats(fulldataset,'place','SI');
sSI_p = cell2mat(stblPlaceSpatialInfo.stable')';
usSI_p = cell2mat(stblPlaceSpatialInfo.unstable')';
grps_pp = [zeros(1,length(sSI_p)), ones(1,length(usSI_p))];
animalColors_pp_stable = nan(length(sSI_p),3);
animalColors_pp_unstable = nan(length(usSI_p),3); 
s = 1;
u = 1;
for a = 1:nAnimals
    animalColors_pp_stable(s:s+stblPlaceSpatialN.stable(a)-1,:) = repmat(colors(a,:),...
        stblPlaceSpatialN.stable(a),1);
    animalColors_pp_unstable(u:u+stblPlaceSpatialN.unstable(a)-1,:) = repmat(colors(a,:),...
        stblPlaceSpatialN.unstable(a),1); 
    
    s = s+stblPlaceSpatialN.stable(a);
    u = u+stblPlaceSpatialN.unstable(a);
end

[stblTimeSpatialInfo,stblTimeSpatialN] = PartitionStats(fulldataset,'time','SI');
sSI_t = cell2mat(stblTimeSpatialInfo.stable')';
usSI_t = cell2mat(stblTimeSpatialInfo.unstable')';
grps_sp = [zeros(1,length(sSI_t)), ones(1,length(usSI_t))];
animalColors_sp_stable = nan(length(sSI_t),3);
animalColors_sp_unstable = nan(length(usSI_t),3); 
s = 1;
u = 1;
for a = 1:nAnimals
    animalColors_sp_stable(s:s+stblTimeSpatialN.stable(a)-1,:) = repmat(colors(a,:),...
        stblTimeSpatialN.stable(a),1);
    animalColors_sp_unstable(u:u+stblTimeSpatialN.unstable(a)-1,:) = repmat(colors(a,:),...
        stblTimeSpatialN.unstable(a),1); 
    
    s = s+stblTimeSpatialN.stable(a);
    u = u+stblTimeSpatialN.unstable(a);
end


%%
fPos = [520 350 300 450];
teal = [0 .5 .5];
purple = [0.5765 0.4392 0.8588];
boxScatterplot([sTI_t,usTI_t],grps_t,'xLabels',{'Stable','Unstable'},...
    'yLabel','Temporal Information [bits]','boxColor',teal,'position',fPos,...
    'circleColors',[animalColors_t_stable;animalColors_t_unstable]);
[~,kp] = kstest2(sTI_t,usTI_t);
tp = ranksum(sTI_t,usTI_t);
title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)]});
 
boxScatterplot([sTI_p,usTI_p],grps_p,'xLabels',{'Stable','Unstable'},...
    'yLabel','Temporal Information [bits]','boxColor',purple,'position',fPos,...
    'circleColors',[animalColors_p_stable;animalColors_p_unstable]);
[~,kp] = kstest2(sTI_p,usTI_p);
tp = ranksum(sTI_p,usTI_p);
title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)]});

boxScatterplot([sSI_p,usSI_p],grps_pp,'xLabels',{'Stable','Unstable'},...
    'yLabel','Spatial Information [bits]','boxColor',purple,'position',fPos,...
    'circleColors',[animalColors_pp_stable;animalColors_pp_unstable]);
[~,kp] = kstest2(sSI_p,usSI_p);
tp = ranksum(sSI_p,usSI_p);
title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)]});

boxScatterplot([sSI_t,usSI_t],grps_sp,'xLabels',{'Stable','Unstable'},...
    'yLabel','Spatial Information [bits]','boxColor',teal,'position',fPos,...
    'circleColors',[animalColors_sp_stable;animalColors_sp_unstable]);
[~,kp] = kstest2(sSI_t,usSI_t);
tp = ranksum(sSI_t,usSI_t);
title({['KS p = ',num2str(kp)], ['T p = ',num2str(tp)]});
