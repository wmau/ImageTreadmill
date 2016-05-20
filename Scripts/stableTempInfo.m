stableI = [];
unstableI = [];
for i=246:254
    [~,a,b] = plotStableTempInfo(MD(247),MD(i),MD(i+1),'time',0);
    stableI = [stableI; a];
    unstableI = [unstableI; b];
end

close all;
figure; 
subplot(1,2,1); hold on;
[~,edges] = histcounts(stableI,linspace(0,max([stableI;unstableI]),20));
[n,b] = hist(unstableI,edges);      %Bin.
n = n./length(unstableI);           %Normalize. 
stairs(b,n,'linewidth',3,'color',[0.7 0.7 0.7]); 
[n,b] = hist(stableI,edges);        %Bin.
n = n./length(stableI);             %Normalize. 
stairs(b,n,'linewidth',3,'color','k');
    xlim([0 max([stableI;unstableI])]);
    set(gca,'ticklength',[0 0]);
    xlabel('Temporal Information [bits/s]');
    ylabel('Proportion'); 
    
subplot(1,2,2); hold on;
[fs,xs] = ecdf(stableI); 
[fn,xn] = ecdf(unstableI); 
plot(xn,fn,'color',[0.7 0.7 0.7],'linewidth',3);
plot(xs,fs,'color','k','linewidth',3);
    xlim([0 max([stableI;unstableI])]);
    legend('Unstable','Stable','location','southeast'); 
    xlabel('Temporal Information [bits/s]'); 
    ylabel('Proportion'); 
    set(gca,'ticklength',[0 0]);
   
[h,p] = kstest2(stableI,unstableI)