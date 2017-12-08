%% Description
% Find the proportion of cells that are either time, place, or time/place
% cells. 
%

%% Code.
clear; 
loadMD;

fulldataset = MD(292:309); 
nSessions = length(fulldataset); 

[pTime,pPlace,pBoth] = deal(nan(nSessions,1)); 
for s=1:nSessions
    pTime(s) = PropCoding(fulldataset(s),'timecells'); 
    pPlace(s) = PropCoding(fulldataset(s),'placecells');
    pBoth(s) = PropCoding(fulldataset(s),'dual'); 
end

figure;
venn([1,mean(pTime),mean(pPlace)],...
    [mean(pTime) mean(pPlace) mean(pBoth) mean(pBoth)],...
    'facecolor',{'k',[0 .5 .5],[.58 .44 .86]});
axis equal;
axis off; 