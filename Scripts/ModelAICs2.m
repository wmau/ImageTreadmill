[~,snks] = find(graphData_p.A);
snks = unique(snks); 

[~,i] = ismember(snks,neurons);
B = size(randcells,1); 

withCellsAIC = zeros(1,length(snks));
noCellsAIC = zeros(1,length(snks)); 
p = zeros(1,length(snks));
randCellsAIC = zeros(B,length(snks));
for n=1:length(snks)
    withCellsAIC(n) = automated{2,i(n)}.ModelCriterion.AIC;
    noCellsAIC(n) = wocells{i(n)}.ModelCriterion.AIC; 
    
    for nn=1:B
        randCellsAIC(nn,n) = randcells{nn,snks(n)}.ModelCriterion.AIC;
    end
    
    p(n) = sum(repmat(withCellsAIC(n),1,B) > randCellsAIC(:,n)')/B;
end

%Jitter for scatter plot. 
wcJitter = 1 - (.1*randn(length(snks),1));
ncJitter = 2 - (.1*randn(length(snks),1)); 
rcJitter = 3 - (.1*randn(length(snks)*B,1)); 

%Groups for boxplotting. 
figure; hold on;
grps = [zeros(1,length(snks)), ones(1,length(snks)), 2*ones(1,length(snks)*B)];
boxplot([withCellsAIC,noCellsAIC,randCellsAIC(:)'],grps,'color','k',...
    'symbol','k','labels',{'Cells','No cells','Random cells'});
scatter([wcJitter; ncJitter; rcJitter],[withCellsAIC'; noCellsAIC'; randCellsAIC(:)],...
    5,'markeredgecolor',[.7 .7 .7]);
set(gca,'ticklength',[0 0]);
ylabel('AIC');

