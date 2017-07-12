%% Description
% Show that trial onsets of time cells (or dual time/place cells) on the
% treadmill and on the track are independent. Do this by correlating "trial
% skewness" of treadmill vs track trial onset. 
%
% Haven't decided whether it's more legit to use time cells or dual
% time/place cells. To switch between the two, go into getAllSkewnesses and
% change the line that reads "neurons =
% AcquireTimePlaceCells(md,'timecells')". 

%% Code.
clear;
loadMD;

fulldataset = MD(292:309);
nSessions = length(fulldataset); 
cellType = 'dual';

[T_Skew,P_Skew] = deal(cell(nSessions,1)); 
for s=1:nSessions
    T_Skew{s} = getAllSkewnesses(fulldataset(s),'cellType',cellType,...
        'rasterType','time');
    P_Skew{s} = getAllSkewnesses(fulldataset(s),'cellType',cellType,....
        'rasterType','place');
end

%Concatenate.
T_Skew_Concat = cell2mat(T_Skew);
P_Skew_Concat = cell2mat(P_Skew);

%Correlation.
[r,p] = corr(T_Skew_Concat,P_Skew_Concat,'rows','complete'); 

figure;
s = scatter(T_Skew_Concat,P_Skew_Concat,50,'k');
set(gca,'tickdir','out','fontsize',12,'linewidth',4);
xlabel('Treadmill within-session trial bias','fontsize',15);
ylabel('Track within-session trial bias','fontsize',15); 
title(['R = ',num2str(r),', P = ',num2str(p)],'fontsize',15);
axis equal;
lsline; 