meanRatio = [];
t = [];
for i=252:255
    cd(MD(i).Location);
    
    load('ProcOut.mat','NumNeurons');
      
    if i~=252
        MakeGraphv4(MD(i));
    end
    
    load('Graphv4.mat','graphData');
    
    if i~=252
        graphData_p = pruneA(MD(i),graphData); 
    end
    
    p = ProgressBar(NumNeurons);
    for neuron=1:NumNeurons
        [ratio,~,~,tmalignedonsets] = SpreadRatio(MD(i),graphData,neuron);
        if length(ratio) > 2
            meanRatio = [meanRatio; mean(ratio)];
            t = [t; mean(cellfun(@mean,tmalignedonsets))];
        else
            meanRatio = [meanRatio; nan];
            t = [t; nan];
        end 
        p.progress;
    end
    
end
p.stop;

figure;
scatter(t,meanRatio,'.');
[r,p] = corr(meanRatio(isfinite(meanRatio)),t(isfinite(meanRatio)),'type','spearman');
title(['R=',num2str(r),', p=',num2str(p)]);