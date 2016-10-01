function [tspikes,ntspikes] = CompareTransOnsetTimes(md)
%
%
%

%% 
    cd(md.Location); 
    load('Pos_align.mat','FT');
    load('TimeCells.mat','TimeCells','T','TodayTreadmillLog');
    load('graphData_p.mat','A'); 
    nNeurons = size(FT,1);
    [~,active] = nNeuronsActiveonTM(md);
    
    %Make treadmill run indices even.
    inds = TodayTreadmillLog.inds; 
    inds = inds(find(TodayTreadmillLog.complete),:);
    inds(:,2) = inds(:,1) + 20*T-1; 
    
    %Get target cells and nontargets. 
    [~,targets] = find(A); 
    targets = unique(targets); 
    nTargets = length(targets); 
    nontargets = setdiff(TimeCells,targets); 
    nNonTargets = length(nontargets);
    
    %Preallocate. 
    tspikes = cell(1,nTargets);
    ntspikes = cell(1,nNonTargets); 
    
    rasters = cell(1,nNeurons);
    for n=active
        rasters{n} = buildRaster(inds,FT,n);
    end
    
    TSPIKES = []; NTSPIKES = [];
    for t=1:nTargets
        triggers = find(A(:,targets(t)))';
        
        tspk = [];
        lspk = [];
        for trig=triggers    
            
            [pRaster,latency,lap] = stripRaster(rasters{trig},rasters{targets(t)});
            
            [spkLap,spkTime] = find(pRaster); 
            [spkLap,order] = sort(spkLap); 
            spkTime = spkTime(order); 
            try
            spkTime = spkTime./20 - latency; 
            catch, keyboard; end
            
            spkalreadythere = ismember(spkTime,tspk);
            lapalreadythere = ismember(spkLap,lspk);
            tspk = [tspk; spkTime(~(spkalreadythere & lapalreadythere))];
            lspk = [lspk; spkLap(~(spkalreadythere & lapalreadythere))];
            
        end
        
        tspikes{t} = tspk;
        TSPIKES = [TSPIKES; tspikes{t}];
    end
    
    for nt=1:nNonTargets
        [~,ntspikes{nt}] = find(rasters{nontargets(nt)});
        ntspikes{nt} = ntspikes{nt}./20;
        NTSPIKES = [NTSPIKES; ntspikes{nt}];
    end
    
    figure;
    [n,b] = hist(TSPIKES,[0:.25:10]); 
    n = n./length(TSPIKES);
    stairs(b,n); 
    hold on;
    [n,b] = hist(NTSPIKES,[0:.25:10]);
    n = n./length(NTSPIKES);
    stairs(b,n); 
    xlabel('Treadmill Time Elapsed [s]'); 
    ylabel('Frequency'); 
    
    figure;
    bar([1,2],[mean(TSPIKES) mean(NTSPIKES)],.4,'facecolor','k'); hold on;
    errorbar(1,mean(TSPIKES),1.96*std(TSPIKES)/sqrt(length(TSPIKES)),'k','linewidth',3)
    errorbar(2,mean(NTSPIKES),1.96*std(NTSPIKES)/sqrt(length(NTSPIKES)),'k','linewidth',3)
    set(gca,'xtick',[1 2],...
        'xticklabel',{'Cell modulated','Treadmill modulated'},...
        'ticklength',[0 0]);
    xlim([0.5 2.5]);
    ylabel('Treadmill Time Elapsed [s]'); 

    keyboard
end