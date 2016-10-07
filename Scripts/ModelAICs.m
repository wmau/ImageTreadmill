[~,sinks] = find(graphData_p.A);
sinks = unique(sinks);

[~,i] = ismember(sinks,neurons);
allAIC = zeros(1,length(i));
nocellsAIC = zeros(1,length(i));
for n=1:length(i)
    allAIC(n) = automated{2,i(n)}.ModelCriterion.AIC;
    %notimeAIC(n) = sNoTime{1,i(n)}.ModelCriterion.AIC;
    nocellsAIC(n) = nocells{3,i(n)}.ModelCriterion.AIC;
end

aJitter = 1-(0.1*randn(length(i),1));
ncJitter = 2-(0.1*randn(length(i),1));

figure; 
subplot(1,2,1); hold on;
% mAll = mean(allAIC);            allSEM = std(allAIC)/sqrt(length(i));
% mNoCells = mean(nocellsAIC);    noCellsSEM = std(nocellsAIC)/sqrt(length(i));
% bar([1,2],[mNoCells,mAll]);
% errorbar(1,mNoCells,noCellsSEM);
% errorbar(2,mAll,allSEM);
grps = [zeros(1,length(i)) ones(1,length(i))];
boxplot([allAIC,nocellsAIC],grps,'color','k','symbol','k','labels',{'Cells','No cells'});
scatter([aJitter; ncJitter],[allAIC'; nocellsAIC'],5,'markeredgecolor',...
    [.7 .7 .7]);
set(gca,'ticklength',[0 0]);
ylabel('AIC');

subplot(1,2,2); 
scatter(allAIC,nocellsAIC,5,'k','filled');
line([-15000 0],[-15000 0],'color','k');
xlabel('Cell Regression AIC');
ylabel('No Cell Regression AIC');

p = signrank(allAIC,nocellsAIC)