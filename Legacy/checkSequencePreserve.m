function checkSequencePreserve(mapMD,base,comp,Ts)
%
%
%

%%
    [normtilemat,sortedPeaks] = msPastalkovaPlot(mapMD,base,comp,Ts);
    close all;
    
    [nTCs,nSessions] = size(sortedPeaks);
    
    diffPeaks = nan(nTCs,nSessions);
    shuffmat = normtilemat; 
    B = 1000;
    shuffdPeaks = nan(nTCs,nSessions,B); 
    for s=2:nSessions
        %Find difference in peaks of the same cell across days. 
        diffPeaks(:,s) = abs(sortedPeaks(:,1) - sortedPeaks(:,s)); 
        
        %Bootstrap. 
        for i=1:B
            shuffmat{s} = normtilemat{s}(randperm(nTCs),:);
            missing = isnan(shuffmat{s}(:,1));
            shuffPeaks = nan(nTCs,1);
            [~,shuffPeaks(~missing,:)] = max(shuffmat{s}(~missing,:),[],2,'includenan');

            shuffPeaks = shuffPeaks./4;
            shuffdPeaks(:,s,i) = abs(sortedPeaks(:,1) - shuffPeaks);
        end
    end
    
    
    null = shuffdPeaks(:);
    figure; hold on;
    [~,edges] = histcounts(null,linspace(0,max(Ts),max(Ts)*4));
    [n,b] = hist(null,edges); 
    n = n./sum(~isnan(null));
    stairs(b,n,'linewidth',3,'color','k'); 
    p = zeros(1,nSessions-1);
    for s=2:nSessions
        [n,b] = hist(diffPeaks(:,s),edges); 
        n = n./sum(~isnan(diffPeaks(:,s))); 
        
        stairs(b,n,'linewidth',3,'color',0.2*s*ones(1,3)); 
        
        [~,p(s-1)] = kstest2(diffPeaks(:,s),null)
    end
    
    xlabel('Absolute Difference from Sequence 1'); 
    ylabel('Proportion'); 
    set(gca,'ticklength',[0 0]);
end