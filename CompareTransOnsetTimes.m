function [ispikes,xspikes] = CompareTransOnsetTimes(md)
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
    ispikes = cell(1,nTargets);
    xspikes = cell(1,nNonTargets); 
    
    rasters = cell(1,nNeurons);
    for n=active
        rasters{n} = buildRaster(inds,FT,n);
    end
    
    iSPIKES = []; xSPIKES = [];
    for t=1:nTargets
        triggers = find(A(:,targets(t)))';
        
        tspk = [];
        lspk = [];
        for trig=triggers    
            
            [pRaster,latency,lap] = stripRaster(rasters{trig},rasters{targets(t)});
            
            [spkLap,spkTime] = find(pRaster); 
            [spkLap,order] = sort(spkLap); 
            spkTime = spkTime(order); 
            spkTime = spkTime./20 -.05 - latency; 
            
            spkalreadythere = ismember(spkTime,tspk);
            lapalreadythere = ismember(spkLap,lspk);
            tspk = [tspk; spkTime(~(spkalreadythere & lapalreadythere))];
            lspk = [lspk; spkLap(~(spkalreadythere & lapalreadythere))];
            
        end
        
        ispikes{t} = tspk;
        iSPIKES = [iSPIKES; ispikes{t}];
    end
    
    for nt=1:nNonTargets
        [~,xspikes{nt}] = find(rasters{nontargets(nt)});
        xspikes{nt} = xspikes{nt}./20 - 0.05;
        xSPIKES = [xSPIKES; xspikes{nt}];
    end
    nXSPIKES = length(xSPIKES); 
    nISPIKES = length(iSPIKES);
    
    figure;
    [n,b] = hist(iSPIKES,[0:.25:10]); 
    n = n./length(iSPIKES);
    stairs(b,n); 
    hold on;
    [n,b] = hist(xSPIKES,[0:.25:10]);
    n = n./length(xSPIKES);
    stairs(b,n); 
    xlabel('Treadmill Time Elapsed [s]'); 
    ylabel('Frequency'); 
    
    %Boxplot stuff.
    figure;
%     bar([1,2],[mean(TSPIKES) mean(NTSPIKES)],.4,'facecolor','k'); hold on;
%     errorbar(1,mean(TSPIKES),1.96*std(TSPIKES)/sqrt(length(TSPIKES)),'k','linewidth',3)
%     errorbar(2,mean(NTSPIKES),1.96*std(NTSPIKES)/sqrt(length(NTSPIKES)),'k','linewidth',3)
%     set(gca,'xtick',[1 2],...
%         'xticklabel',{'Cell modulated','Treadmill modulated'},...
%         'ticklength',[0 0]);
%     xlim([0.5 2.5]);
    xJitter = 1-(0.1*randn(nXSPIKES,1));
    iJitter = 2-(0.1*randn(nISPIKES,1));
    grps = [zeros(length(xSPIKES),1);ones(length(iSPIKES),1)];

    hold on;
    scatter([xJitter; iJitter],[xSPIKES; iSPIKES],3,...
        'markeredgecolor',[.7 .7 .7]);
    boxplot([xSPIKES;iSPIKES],grps,'color','k','whisker',0,'symbol','k',...
        'Labels',{'Treadmill','Cell'});
    ylabel('Treadmill Time Elapsed [s]'); 
    xlabel('Modulation');
    set(gca,'ticklength',[0 0]);
end