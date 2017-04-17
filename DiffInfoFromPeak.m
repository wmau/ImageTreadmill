function [alignedDiffs,stable,unstable,dayDelta,alignedDays] = DiffInfoFromPeak(mds,cellType,infoType,varargin)
%[alignedDiffs,stable,unstable] = DiffInfoFromPeak(mds,cellType,infoType,varargin)
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
    end
    DATA = CompileMultiSessionData(mds,{cellTypeString});
    
    %Map them.
    codingCells = DATA.(cellTypeString);
    map = msMatchMultiSessionCells(mds,codingCells);
    
    %Get total number of neurons. 
    nNeurons = size(map,1);
    
%% Get all informations.
    [stability,I] = deal(nan(size(map)));
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
        
        %Dump into matrix.
        I(mapped,s) = MI(neurons(mapped));
        
        %Determine stability.
        switch cellType
            case 'time', nCodingCells = length(getTimeCells(mds(s))); 
            case 'place',nCodingCells = length(getPlaceCells(mds(s),.01)); 
        end       
        crit = .01/nCodingCells;
        
        %Perform cross-days correlation.
        if s~=nSessions
            switch cellType
                case 'time', corrStats = CorrTrdmllTrace(mds(s),mds(s+1),neurons(mapped));
                case 'place',corrStats = CorrPlaceFields(mds(s),mds(s+1),neurons(mapped));
            end  
            
            stability(mapped,s) = corrStats(neurons(mapped),2) < crit;
        end
    end

%% Align to max.
    %Identify the max information and day on which it occurred.
    [peakI,peakDay] = max(I,[],2); 
    diffFromPeak = I - peakI;
    propOfPeak = I./peakI;
    
    %Max day difference.
    maxDayDelta = nSessions-1; %Plus or minus.
    dayDelta = [-maxDayDelta:maxDayDelta];
    nDays = length(dayDelta);
    
    %Make matrix of aligned information.
    [alignedDays,alignedDiffs] = deal(nan(nNeurons,nDays));
    for s=1:nSessions
        %Distance to peak day.
        d = s-peakDay; 
            
        %Get column where information difference goes. 
        [~,alignedColInd] = ismember(d,dayDelta);
        
        %Dump information differentials into matrix. 
        for n=1:nNeurons
            alignedDiffs(n,alignedColInd(n)) = propOfPeak(n,s);  
            alignedDays(n,alignedColInd(n)) = d(n);
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