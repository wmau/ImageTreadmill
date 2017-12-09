clear
loadMD;
fulldataset = MD(292:309);

saveBool = false;
folder = 'C:\Users\William Mau\Documents\Projects\Time Cell Imaging Summer 2015 -\Paper\Figures\Supplementals\Alternative Info Stability Analyses';
TimeTI = fullfile(folder,'Correlate Rho to Info_Time Stability given TI');
TimeSI = fullfile(folder,'Correlate Rho to Info_Time Stability given SI');
TimeFR = fullfile(folder,'Correlate Rho to Info_Time Stability given FR');
PlaceSI = fullfile(folder,'Correlate Rho to Info_Place Stability given SI');
PlaceTI = fullfile(folder,'Correlate Rho to Info_Place Stability given TI');
PlaceFR = fullfile(folder,'Correlate Rho to Info_Place Stability given FR'); 
TimeHistogram = fullfile(folder,'Temporal Stability Rhos');
PlaceHistogram = fullfile(folder,'Spatial Stability Rhos');

if saveBool
    c = input('Saving set to true. Are you sure you want to continue? (y/n)','s');

    if ~strcmp(c,'y')
        saveBool = false;
    end
end

%%
CORRS = corrInfoStability(fulldataset,'time','ti');
if saveBool, print(TimeTI,'-dpdf'); end

% corrInfoStability(fulldataset,'time','si');
% if saveBool, print(TimeSI,'-dpdf'); end
% 
% corrInfoStability(fulldataset,'time','fr');
% if saveBool, print(TimeFR,'-dpdf'); end

% figure;
% histogram(CORRS,50,'normalization','probability','edgecolor','none');
% set(gca,'tickdir','out','linewidth',3);
% xlabel('Temporal Stability Rho');
% ylabel('Proportion');
% if saveBool, print(TimeHistogram,'-dpdf'); end

%%
CORRS = corrInfoStability(fulldataset,'place','si');
if saveBool, print(PlaceSI,'-dpdf'); end

% corrInfoStability(fulldataset,'place','ti');
% if saveBool, print(PlaceTI,'-dpdf'); end
% 
% corrInfoStability(fulldataset,'place','fr');
% if saveBool, print(PlaceFR,'-dpdf'); end
% 
% figure;
% histogram(CORRS,50,'normalization','probability','edgecolor','none');
% set(gca,'tickdir','out','linewidth',3);
% xlabel('Spatial Stability Rho');
% ylabel('Proportion');
% if saveBool, print(PlaceHistogram,'-dpdf'); end
