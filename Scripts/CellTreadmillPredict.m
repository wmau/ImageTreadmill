meanRatio = [];
medRatio = [];
pRatio = [];
t = [];
for i=215%:2:221
    cd(MD(i).Location);
    
    load('FinalOutput.mat','FT');
    NumNeurons = size(FT,1);
    load('Graphv4.mat','graphData');
     
    p = ProgressBar(NumNeurons);
    for neuron=1:NumNeurons
        el = [1:NumNeurons]; el(neuron) = [];
        [ratio,~,~,tmalignedonsets] = SpreadRatio(MD(i),graphData,neuron);
            meanRatio = [meanRatio; nanmean(ratio)];
            medRatio = [medRatio; nanmedian(ratio)];
            pRatio = [pRatio; sum(ratio < 1)/length(ratio)];
            t = [t; nanmean(cellfun(@mean,tmalignedonsets))];

        p.progress;
    end
    
    nTrigger = sum(graphData_p.A);
    pTrigger = nTrigger./sum(~isnan(graphData_p.prune_p));
    pTrigger = pTrigger';
end
p.stop;

figure;
subplot(4,1,1);
scatter(t,meanRatio,'.');
[r,p] = corr(meanRatio(isfinite(meanRatio)),t(isfinite(meanRatio)),'type','spearman');
title(['R=',num2str(r),', p=',num2str(p)]); ylabel('Mean Ratio'); 

subplot(4,1,2);
scatter(t,medRatio,'.');
[r,p] = corr(medRatio(isfinite(medRatio)),t(isfinite(medRatio)),'type','spearman');
title(['R=',num2str(r),', p=',num2str(p)]); ylabel('Median Ratio');

subplot(4,1,3);
scatter(t,pRatio,'.');
[r,p] = corr(pRatio(isfinite(pRatio)),t(isfinite(pRatio)),'type','spearman');
title(['R=',num2str(r),', p=',num2str(p)]); ylabel('P(Ratio < 1)');

subplot(4,1,4);
scatter(t,pTrigger);
[r,p] = corr(pTrigger(isfinite(pTrigger)),t(isfinite(pTrigger)),'type','spearman');
title(['R=',num2str(r),', p=',num2str(p)]); ylabel('P(Internally modulated)');
