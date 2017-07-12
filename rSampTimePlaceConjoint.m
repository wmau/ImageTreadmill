function [dual_pct,mean_rSampDual_pct] = rSampTimePlaceConjoint(md)
%
%
%

%%
    cd(md.Location); 
    B = 1000;
    
    %Get time and place cells, as well as conjunctive time/place cells. 
    TCs = AcquireTimePlaceCells(md,'timecells');
    PCs = AcquireTimePlaceCells(md,'placecells');
    dual = AcquireTimePlaceCells(md,'dual');
    
    %Get number of time and place cells. 
    nTCs = length(TCs);
    nPCs = length(PCs); 
    
    load('FinalOutput.mat','NumNeurons');
    N = length(union(TCs,PCs));
    
    rSampDual_pct = nan(B,1);
    for i=1:B
        rSampTime = randsample(1:NumNeurons,nTCs);
        rSampPlace = randsample(1:NumNeurons,nPCs); 
        rSampDual = intersect(rSampTime,rSampPlace);
        
        rSampDual_pct(i) = length(rSampDual)/N; 
    end
    
    dual_pct = length(dual)/N;
    p = sum(dual_pct > rSampDual_pct)/B;
   
    mean_rSampDual_pct = mean(rSampDual_pct); 
%     figure;
%     histogram(rSampDual_pct,'edgecolor','none');
%     yLims = get(gca,'ylim');
%     line([dual_pct dual_pct],[0 yLims(2)],'color','r'); 
%     set(gca,'tickdir','out','linewidth',4,'fontsize',12);
%     xlabel('Proportion conjoint time/place cells'); 
%     ylabel('Count'); 
%     title(['P = ',num2str(p)]);
end
    
    