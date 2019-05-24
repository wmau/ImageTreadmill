function plotTraceWithTreadmillEpochs(md,varargin)
%
%
%

%%
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addParameter('timecells',[]);
    
    p.parse(md,varargin{:});
    timecells = p.Results.timecells;
    if isempty(timecells)
        timecells = getTimeCells(md); 
    end
    nTimeCells = length(timecells); 
%%
    load(fullfile(md.Location,'Pos_align.mat'),'RawTrace');
    traces = RawTrace;
    load(fullfile(md.Location,'TimeCells.mat'),'TodayTreadmillLog','T'); 
    t = (1:size(traces,2))/20; 
    
    onTreadmill = TrimTrdmllInds(TodayTreadmillLog,T);
    nTrials = size(onTreadmill,1); 
    
    keepgoing = true;
    i = 1;
    while keepgoing
        f = figure;
        
        plot(t,traces(timecells(i),:));
        hold on;
        
        for thisTrial = 1:nTrials
            epoch = [onTreadmill(thisTrial,1):onTreadmill(thisTrial,2)];
            plot(t(epoch),traces(timecells(i),epoch),'r');
        end
        
        title(['Neuron #',num2str(timecells(i))]);
        xlabel('Time (s)');
        ylabel('DF/F');
        make_plot_pretty(gca); 
        
        xlim([t(1) t(end)]);
        
        [keepgoing,i] = scroll(i,nTimeCells,f);
        close all;
    end
end