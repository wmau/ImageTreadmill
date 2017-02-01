clear
loadMD; 
fulldataset = MD(292:309);

teal = [0 .5 .5];
purple = [0.5765 0.4392 0.8588];

folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals';
time = fullfile(folder,'Info vs Days Stable (Time Stability)');
place = fullfile(folder,'Info vs Days Stable (Place Stability)');

saveBool = true;
if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');

    if ~strcmp(c,'y')
        saveBool = false;
    end
end

figure;
InformationOverDaysStable(fulldataset,'time','ti'); 
hold on;
InformationOverDaysStable(fulldataset,'time','si');
e = get(gca,'Children');
e(1).XData = e(1).XData + .2;
set(gca,'xtick',[0:4]);
e(1).Color = purple;
e(2).Color = teal;
legend({'Temporal','Spatial'},'location','northwest');
if saveBool
    print(time,'-dpdf');
end

figure;
InformationOverDaysStable(fulldataset,'place','si'); 
hold on;
InformationOverDaysStable(fulldataset,'place','ti');
e = get(gca,'Children');
e(1).XData = e(1).XData + .2;
set(gca,'xtick',[0:4]);
e(1).Color = teal;
e(2).Color = purple;
legend({'Spatial','Temporal'},'location','northwest');
if saveBool
    print(place,'-dpdf');
end