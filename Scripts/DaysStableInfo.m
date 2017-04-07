clear
loadMD; 
fulldataset = MD(292:309);

teal = [0 .5 .5];
purple = [0.5765 0.4392 0.8588];

folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals';
time = fullfile(folder,'Info vs Days Stable (Time Stability)');
place = fullfile(folder,'Info vs Days Stable (Place Stability)');

saveBool = false;
if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');

    if ~strcmp(c,'y')
        saveBool = false;
    end
end
%%
figure('Position',[520   380   290   420]); hold on
[ttiDuration,ttiStats] = InformationOverDaysStable(fulldataset,'time','ti',teal);    
[tsiDuration,tsiStats] = InformationOverDaysStable(fulldataset,'time','si',purple);

e = get(gca,'Children');
%e(1).XData = e(1).XData + .2;
set(gca,'xtick',[0:1:4],'linewidth',4,'tickdir','out','fontsize',15);
e(1).Color = purple;
e(2).Color = teal;
legend({'Temporal','Spatial'},'location','northwest');
if saveBool
    print(time,'-dpdf');
end

figure;
[tbl,stats,comps] = ANOVAit(ttiDuration,ttiStats,tsiDuration,tsiStats);

%%
figure('Position',[520   380   290   420]); hold on;
[ssiDuration,ssiStats] = InformationOverDaysStable(fulldataset,'place','si',purple); 
[stiDuration,stiStats] = InformationOverDaysStable(fulldataset,'place','ti',teal);

e = get(gca,'Children');
%e(1).XData = e(1).XData + .2;
set(gca,'xtick',[0:1:4],'linewidth',4,'tickdir','out','fontsize',15);
e(1).Color = teal;
e(2).Color = purple;
legend({'Spatial','Temporal'},'location','northwest');
if saveBool
    print(place,'-dpdf');
end

figure;
[tbl,stats,comps] = ANOVAit(ssiDuration,ssiStats,stiDuration,stiStats);
    
function [tbl,stats,comps] = ANOVAit(sameDuration,sameStats,diffDuration,diffStats)
%
%
%

%% 
    dimension = [ones(length(sameStats),1); zeros(length(diffStats),1)];
    duration = [sameDuration; diffDuration];
    X = [sameStats; diffStats]; 
    grps = {dimension,duration};
    [~,tbl,stats] = anovan(X,grps,'model','full','varnames',{'Dimension','Duration'},'display','on');
    comps = multcompare(stats,'dimension',[1,2],'display','on');
end