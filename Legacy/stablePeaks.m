stablePeak = [];
unstablePeak = [];
crit = 'time';
for i=243:246    
    load(fullfile(MD(i).Location,'TimeCells.mat'),'TimeCells');
    [pCorr,tCorr] = PlaceTimeCorr(MD(244),MD(i),MD(i+1),TimeCells);
    
    if strcmp(crit,'time')
        sig = tCorr(:,2) < 0.05; 
        nsig = tCorr(:,2) > 0.05; 
    else 
        sig = pCorr(:,2) < 0.05;
        nsig = pCorr(:,2) > 0.05;
    end
    
    t = getTimePeak(MD(i)); 
    stablePeak = [stablePeak; t(sig)];
    unstablePeak = [unstablePeak; t(nsig)];
end

figure;
subplot(1,2,1); hold on;
[~,edges] = histcounts(stablePeak,linspace(0,max([stablePeak;unstablePeak]),20));
[n,b] = hist(unstablePeak,edges); 
n = n./length(unstablePeak); 
stairs(b,n,'linewidth',3,'color',[0.7 0.7 0.7]); 
[n,b] = hist(stablePeak,edges);        %Bin.
n = n./length(stablePeak);             %Normalize. 
stairs(b,n,'linewidth',3,'color','k');
    xlim([0 max([stablePeak;unstablePeak])]);
    set(gca,'ticklength',[0 0]);
    xlabel('Peak Tuning [s]');
    ylabel('Proportion'); 
    
subplot(1,2,2); hold on;
[fs,xs] = ecdf(stablePeak); 
[fn,xn] = ecdf(unstablePeak); 
plot(xn,fn,'color',[0.7 0.7 0.7],'linewidth',3);
plot(xs,fs,'color','k','linewidth',3);
    xlim([0 max([stablePeak;unstablePeak])]);
    legend('Unstable','Stable','location','southeast'); 
    xlabel('Peak Tuning [s]'); 
    ylabel('Proportion'); 
    set(gca,'ticklength',[0 0]);

[h,p] = kstest2(stablePeak,unstablePeak)