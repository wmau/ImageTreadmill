animal = 'GCamp6f_45_treadmill';
%date = '11_20_2015';
%RecordStartTime = '10:47:51 AM';

date = '11_19_2015';
RecordStartTime = '4:42:32 PM';

TodayTreadmillLog = getTodayTreadmillLog(animal,date);
TodayTreadmillLog = AlignTreadmilltoTracking(TodayTreadmillLog,RecordStartTime);

load('ProcOut.mat','FT');

[ratebylap,bin] = getLapResponses(animal,date,2,FT,TodayTreadmillLog);

activelaps = nan(size(ratebylap,3),1);
for i=1:size(ratebylap,3)
    activelaps(i) = sum(any(ratebylap(:,:,i),2));
end

%Define number of laps criterion here.
crit = 15; 
good = find(activelaps > crit); 

for i=2:2:length(good)
    figure;
    img = imagesc([0:2:max(bin)],...
        [1:5:sum(TodayTreadmillLog.complete)],...
        ratebylap(:,:,good(i)));
    set(img,'AlphaData',~isnan(ratebylap(:,:,good(i))));
    colormap gray;
    xlabel('Time [s]');
    ylabel('Laps');
end