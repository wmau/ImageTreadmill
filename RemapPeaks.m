function RemapPeaks(mapMD,base,comp,Ts)
%
%
%

%% 
    tResolution = 0.25; 
    nBins = Ts(1)/tResolution;
    nBins2 = Ts(2)/tResolution;
    t = linspace(0,Ts(1),nBins)';
    t2 = linspace(0,Ts(2),nBins2)';
    
    [sametuning,TIMECELLS,CURVES,MAP,MAPcols] = TimeCellRemapRate(mapMD,base,comp,Ts);
    
    stable = TIMECELLS{1}(sametuning(:,2) == 1); 
    unstable = TIMECELLS{1}(sametuning(:,2) ~= 1); 
    
    [~,sbinpks] = cellfun(@max,CURVES{1}.tuning(stable));
    [~,ubinpks] = cellfun(@max,CURVES{1}.tuning(unstable));
    
    stable_tpks = t(sbinpks);   
    unstable_tpks = t(ubinpks);
    
    [sN,sBins] = hist(stable_tpks,[0:Ts(1)]);   
    [uN,uBins] = hist(unstable_tpks,[0:Ts(1)]); 

    [~,p] = kstest2(stable_tpks,unstable_tpks); 
    
    sN = sN./sum(sN); 
    uN = uN./sum(uN);

    figure;
    stairs(sBins,sN,'linewidth',4,'color','k');             hold on;
    stairs(uBins,uN,'linewidth',4,'color',[0.8 0.8 0.8]);   hold off; 
    set(gca,'xtick',0:2:Ts(1));
    set(gca,'ticklength',[0 0]);
    xlabel('Time [s]'); 
    ylabel('Proportion'); 
    legend({'Stable','Unstable'}); 
    title(['p=',num2str(p)]);
    set(gca,'fontsize',15); 

    stable_ind = MAP(ismember(MAP(:,MAPcols(1)),stable),MAPcols(2));    
    unstable_ind = MAP(ismember(MAP(:,MAPcols(1)),unstable),MAPcols(2));
    
    stable_good = stable_ind ~= 0; 
    unstable_good = unstable_ind ~= 0; 
    
    stable_ind(stable_ind==0) = []; unstable_ind(unstable_ind==0) = []; 
    
    [~,sbinpks2] = cellfun(@max,CURVES{2}.tuning(stable_ind));
    [~,ubinpks2] = cellfun(@max,CURVES{2}.tuning(unstable_ind));
    
    stable_tpks2 = t2(sbinpks2); 
    unstable_tpks2 = t2(ubinpks2); 
    
    figure;
    scatter(stable_tpks(stable_good),stable_tpks2,2000,'k.');                       hold on;
    scatter(unstable_tpks(unstable_good),unstable_tpks2,2000,[0.8 0.8 0.8],'.');    hold off;
    xlabel('Peak on Day 1 [s]'); ylabel('Peak on Day 2 [s]'); 
    legend({'Stable','Unstable'},'location','best'); 
    set(gca,'fontsize',15); 
    set(gca,'ticklength',[0 0]);
    

end