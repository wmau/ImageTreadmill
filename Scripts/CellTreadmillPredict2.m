[~,target] = find(graphData.A); 
target = unique(target)'; 

target = 1:1339;

t = [];
c = 1;
pCellBetter = zeros(length(target),1);
p = ProgressBar(length(target));
for neuron = target
    el = target;
    el(el==neuron) = [];
    [ratio,~,~,tmalignedonsets] = SpreadRatio(MD(252),graphData,neuron,'edgelist',el);
    
    pCellBetter(c) = sum(ratio < 1)/length(ratio);
    t = [t; nanmean(cellfun(@nanmean, tmalignedonsets))];
    
    c=c+1;
    p.progress;
end
p.stop;

figure;
scatter(t,pCellBetter,'.');
[r,p] = corr(pCellBetter(isfinite(pCellBetter)),t(isfinite(pCellBetter)),'type','spearman','rows','complete');
title(['R=',num2str(r),', p=',num2str(p)]);