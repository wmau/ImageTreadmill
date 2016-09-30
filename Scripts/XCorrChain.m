path = AllPaths{7}{l10(15)};
nNeurons = length(path); 
colors = rand(nNeurons,3);
md = MD(215);

loadMD;
cd(md.Location);

%Load and normalize traces.
load('Pos_align.mat','rawtrace');
FT = rawtrace; clear rawtrace; 
[~,nFrames] = size(FT);
FT = FT./repmat(max(FT,[],2),1,nFrames);

load('TimeCells.mat','T','TodayTreadmillLog');

%Make treadmill run indices even.
inds = TodayTreadmillLog.inds; 
inds = inds(find(TodayTreadmillLog.complete),:);
inds(:,2) = inds(:,1) + 20*T-1; 

rasters = cell(1,nNeurons);
for n=1:nNeurons
    rasters{n} = buildRasterTrace(inds,FT,path(n));
end

nLaps = size(rasters{1},1);
for trig=1:nNeurons
    thisTrig = path(trig);
    
    for targ=1:nNeurons
        thisTarg = path(targ); 
        
        %XCorr.
        [R.unsmoothed{thisTrig,thisTarg},lags] = xcorr_by_laps(rasters{trig},rasters{targ});
        R.smoothed{thisTrig,thisTarg} = mean(R.unsmoothed{thisTrig,thisTarg});
        
        %Confidence interval.                        
        R.CI{thisTrig,thisTarg} = 1.96*std(R.unsmoothed{thisTrig,thisTarg})/sqrt(nLaps);
    end
end

lags = lags./20;

figure('Position',[390 70 780 730]); 
line1.col = {'k'};
line2.col = {[.6 .6 .6]};
for i=nNeurons:-1:2
    subplot(3,3,i-1);
    mseb(lags,R.smoothed{path(i),path(i-1)},R.CI{path(i),path(i-1)},line1,1);
    
    if i>2
        mseb(lags,R.smoothed{path(i),path(i-2)},R.CI{path(i),path(i-2)},line2,1);
    end
    
    if i==5
        ylabel('Cross-Correlation Value');
    end
    
    if i==9
        xlabel('Time [s]');
    end
    
    ylim([-.3 .65]);
    title(['\color[rgb]{',num2str(colors(i-1,:)),'}',num2str(path(i-1)),...
        ' \color{black}x ',...
        '\color[rgb]{',num2str(colors(i,:)),'}',num2str(path(i))]);
end

figure;
PlotPath(md,path,colors);