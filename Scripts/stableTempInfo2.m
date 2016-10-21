stableI = [];
unstableI = [];
for i=215:2:219
    [~,a,b] = plotStableTempInfo(MD(215),MD(i),MD(i+2),'time',0);
    stableI = [stableI; a];
    unstableI = [unstableI; b];
end

%Means
mStableI = mean(stableI);
mUnstableI = mean(unstableI);

%Lengths
nStable = length(stableI);
nUnstable = length(unstableI);

%SEM
N = 1000;
bootstat = sort(bootstrp(N,@mean,stableI));
semStableI = [bootstat(round(N*.95))-mStableI, bootstat(round(N*.05))-mStableI];
bootstat = sort(bootstrp(N,@mean,unstableI));
semUnstableI = [bootstat(round(N*.95))-mUnstableI, bootstat(round(N*.05))-mUnstableI];

n = 3;

%Jitter
sJitter = 1 - (0.2 - 0.4*rand(nStable,1));
uJitter = n - (0.2 - 0.4*rand(nUnstable,1));

close all;
figure; 
subplot(1,2,1);
bar([1;n],[mStableI; mUnstableI],'facecolor','w','edgecolor','k','linewidth',2);
hold on;
scatter([sJitter; uJitter],[stableI; unstableI],15,'o',...
    'markeredgealpha',.4,'markeredgecolor',[0.7 0.7 0.7]);
errorbar(1,mStableI,semStableI(1),semStableI(2),'k','linewidth',2);
errorbar(n,mUnstableI,semUnstableI(1),semUnstableI(2),'k','linewidth',2);
set(gca,'xtick',[1 n],'xticklabel',{'Stable','Unstable'},...
    'linewidth',3,'ticklength',[0 0],'fontsize',16);
ylabel('Temporal Information [bits/s]');

[~,p] = kstest2(stableI,unstableI)

stableI = [];
unstableI = [];
for i=215:2:219
    [~,a,b] = plotStableTempInfo(MD(215),MD(i),MD(i+2),'place',0);
    stableI = [stableI; a];
    unstableI = [unstableI; b];
end

%Means
mStableI = median(stableI);
mUnstableI = median(unstableI);

%Lengths;
nStable = length(stableI);
nUnstable = length(unstableI);

%SEM
semStableI = std(stableI)/sqrt(nStable);
semUnstableI = std(unstableI)/sqrt(nUnstable);

sJitter = 1 - (0.05 - 0.1*rand(nStable,1));
uJitter = n - (0.05 - 0.1*rand(nUnstable,1));

subplot(1,2,2);
bar([1;n],[mStableI; mUnstableI],'facecolor','w','edgecolor','k','linewidth',2);
hold on;
hs = scatter([sJitter; uJitter],[stableI; unstableI],15,'o','markeredgecolor',[0.7 0.7 0.7]);
errorbar(1,mStableI,semStableI,'k','linewidth',2);
errorbar(n,mUnstableI,semUnstableI,'k','linewidth',2);
set(gca,'xtick',[1 n],'xticklabel',{'Stable','Unstable'},...
    'linewidth',3,'ticklength',[0 0],'fontsize',16);
[~,p] = kstest2(stableI,unstableI)
