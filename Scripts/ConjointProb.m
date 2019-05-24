%% Description
% Show that there is no difference in the number of conjoint time/place
% cells than if you randomly sample time and place cells from chance. 
%

%% Code.
clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];   
nSessions = length(fulldataset); 

[empPct,randPct] = deal(nan(nSessions,1)); 
for s=1:nSessions
    [empPct(s),randPct(s)] = rSampTimePlaceConjoint(fulldataset(s)); 
end

figure('Position',[680 620 240 355]);
hold on;
bar([1,2],[mean(empPct),mean(randPct)],'facecolor',[0 .5 .5],'linewidth',4);
errorbar(1,mean(empPct),std(empPct)/sqrt(nSessions),'color','k',...
    'linewidth',4);
errorbar(2,mean(randPct),std(randPct)/sqrt(nSessions),'color','k',...
    'linewidth',4);
for s=1:nSessions
    plot([1,2],[empPct(s),randPct(s)],'color',[.7 .7 .7]);
end
set(gca,'linewidth',4,'tickdir','out','linewidth',4,'xtick',[1 2],'xticklabel',...
    {'Conjoint','Random'},'fontsize',12);
xlim([.5 2.5]);
ylabel('Proportion');

p = ranksum(randPct,empPct)