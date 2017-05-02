function [alignedDiffs,stable,unstable,dayDelta,alignedDays,alignedOtherI] = DiffInfoFromPeak(mds,cellType,infoType,varargin)
%[alignedDiffs,stable,unstable,dayDelta,alignedDays,alignedOtherI] = DiffInfoFromPeak(mds,cellType,infoType,varargin)
%
%

%%
    p = inputParser;
    p.addRequired('mds',@(x) isstruct(x)); 
    p.addParameter('plotit',true,@(x) islogical(x));
    
    p.parse(mds,varargin{:});
    plotit = p.Results.plotit;
    
%% Get all coding cells.
    nSessions = length(mds); 
    switch cellType
        case 'time', cellTypeString = 'timecells'; c = [0 .5 .5];
        case 'place',cellTypeString = 'placecells'; c = [.58 .44 .86];
        case 'dual',cellTypeString = 'dualcells'; c = [1 1 1];
    end
    DATA = CompileMultiSessionData(mds,{cellTypeString});
    
    %Map them.
    codingCells = DATA.(cellTypeString);
    map = msMatchMultiSessionCells(mds,codingCells);
    
    %Get total number of neurons. 
    nNeurons = size(map,1);
    
%% Get all informations.
    [stability,I,otherI] = deal(nan(size(map)));
    for s=1:nSessions
        cd(mds(s).Location);
        
        %Load information.
        switch infoType
            case 'ti', load('TemporalInfo.mat','MI'); 
            case 'si', load('SpatialInfo.mat','MI');
        end
        
        %Get neurons this session whose TI we're extracting. 
        neurons = map(:,s); 
        mapped = find(neurons>0);
        goodCells = neurons(mapped);
        
        %Dump into matrix.
        I(mapped,s) = MI(goodCells);
        
        %Get information from other dimension.
        switch infoType
            case 'ti', load('SpatialInfo.mat','MI');
            case 'si', load('TemporalInfo.mat','MI');
        end
        otherI(mapped,s) = MI(goodCells); 
        
        %Perform cross-days correlation to determine stability.
        if s~=nSessions
            switch cellType
                case 'time', codingCells = getTimeCells(mds(s)); 
                case 'place',codingCells = getPlaceCells(mds(s),.01); 
                case 'dual', codingCells = AcquireTimePlaceCells(mds(s),'dual');
            end       
            codingCells = EliminateUncertainMatches([mds(s),mds(s+1)],codingCells);
            nCodingCells = length(codingCells);
            crit = .01/nCodingCells;
        
            existsAcrossOneDay = EliminateUncertainMatches([mds(s),mds(s+1)],neurons(mapped));
            switch cellType
                case 'time', corrStats = CorrTrdmllTrace(mds(s),mds(s+1),goodCells);
                case 'place',corrStats = CorrPlaceFields(mds(s),mds(s+1),goodCells);
                case 'dual', corrStats = CorrTrdmllTrace(mds(s),mds(s+1),goodCells);
            end  
            
            [~,exists] = ismember(neurons,existsAcrossOneDay);
            exists(exists==0) = [];
            stability(exists,s) = 0;
            stability(mapped,s) = corrStats(goodCells,2) < crit;
        end
    end

%% Align to max.
    %Identify the max information and day on which it occurred.
    [peakI,peakDay] = max(I,[],2); 
    otherIPeak = max(otherI,[],2);          %Peaks for other dimension's information.
    
    %Get proportions of the peak.
    propofOtherIPeak = otherI./otherIPeak;
    propOfPeak = I./peakI;
    
%     propofOtherIPeak = zscore(otherI,[],2);
%     propOfPeak = zscore(I,[],2);
     
    %Or percent change.
%     pctChange = I(:,2:end)./I(:,1:end-1);
%     pctChange = pctChange - 1;
%     pctChange = [nan(size(pctChange,1),1), pctChange];
%     
%     pctOtherChange = otherI(:,2:end)./otherI(:,1:end-1);
%     pctOtherChange = pctOtherChange - 1;
%     pctOtherChange = [nan(size(pctOtherChange,1),1), pctOtherChange];
    
    %Max day difference.
    maxDayDelta = nSessions-1; %Plus or minus.
    dayDelta = [-maxDayDelta:maxDayDelta];
    nDays = length(dayDelta);
    
    %Make matrix of aligned information.
    [alignedDays,alignedDiffs,alignedOtherI] = deal(nan(nNeurons,nDays));
    for s=1:nSessions
        %Distance to peak day.
        d = s-peakDay; 
            
        %Get column where information difference goes. 
        [~,alignedColInd] = ismember(d,dayDelta);
        
        %Dump information differentials into matrix. 
        for n=1:nNeurons
            alignedDiffs(n,alignedColInd(n)) = propOfPeak(n,s);  
            alignedDays(n,alignedColInd(n)) = d(n);
            
            alignedOtherI(n,alignedColInd(n)) = propofOtherIPeak(n,s);
        end     
        
    end
    
    %Determine whether cell was stable for day after the peak.
    stableAfterPeak = nan(nNeurons,1);
    for n=1:nNeurons
        try 
            stableAfterPeak(n) = stability(n,peakDay(n));
        catch 
        end
    end

    %Separate stable vs unstable cells.
    stable = find(stableAfterPeak>0);
    unstable = find(stableAfterPeak==0);
    
    if plotit
        %Get differentials.
        stableDiffs = alignedDiffs(stable,:);
        unstableDiffs = alignedDiffs(unstable,:);

        %Get means.
        meanDiffs = nanmean(alignedDiffs);
        SEM = nanstd(alignedDiffs)./sqrt(nNeurons);

        %Plot.
        figure;
        errorbar(meanDiffs,SEM,'linewidth',4,'color','k'); hold on;
        errorbar(nanmean(stableDiffs),nanstd(stableDiffs)./sqrt(length(stable)),...
            'linewidth',4,'color',c);
        errorbar(nanmean(unstableDiffs),nanstd(unstableDiffs)./sqrt(length(unstable)),...
            'linewidth',4','color','r');
        set(gca,'linewidth',4,'fontsize',12,'tickdir','out','xticklabel',dayDelta);
        xlabel('Days from Peak','fontsize',15);
        ylabel('Proportion of Info.','fontsize',15);

    end
end