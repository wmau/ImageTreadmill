[~,sinks] = find(graphData_p.A);
sinks = unique(sinks);

[~,i] = ismember(sinks,neurons);
for n=1:length(i)
    allAIC(n) = sAll{1,i(n)}.ModelCriterion.AIC;
    notimeAIC(n) = sNoTime{1,i(n)}.ModelCriterion.AIC;
    nocellsAIC(n) = sNoCells{1,i(n)}.ModelCriterion.AIC;
end

figure; 
subplot(1,2,1); hold on;
mAll = mean(allAIC);            allSEM = std(allAIC)/sqrt(length(i));
mNoCells = mean(nocellsAIC);    noCellsSEM = std(nocellsAIC)/sqrt(length(i));
bar([1,2],[mNoCells,mAll]);
errorbar(1,mNoCells,noCellsSEM);
errorbar(2,mAll,allSEM);

subplot(1,2,2); 
scatter(allAIC,nocellsAIC,'.');
refline;
p = signrank(allAIC,nocellsAIC)