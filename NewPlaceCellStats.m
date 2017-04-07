function [newPCStatS1,r1,stat1,newPCStatS2,r2,stat2] = NewPlaceCellStats(base,comp,statType,varargin)
%
%
%

%% Parse inputs. 
    p = inputParser; 
    p.addRequired('base',@(x) isstruct(base));
    p.addRequired('comp',@(x) isstruct(comp)); 
    p.addRequired('statType',@(x) ischar(x));
    p.addParameter('plotit',true,@(x) islogical(x));
    
    p.parse(base,comp,statType,varargin{:});
    
    plotit = p.Results.plotit;

%%
    
    [newPCsS1,newPCsS2] = getNewPlaceCells(base,comp);
    S1PCs = getPlaceCells(base,.01);
    S2PCs = getPlaceCells(comp,.01);
    S1TCs = getTimeCells(base);
    S2TCs = getTimeCells(comp);
    cd(base.Location);

    binsize = .1;
    switch statType
        case 'fr'
            load('Pos_align.mat','PSAbool');
            [n,f] = size(PSAbool);
            d = diff([zeros(n,1) PSAbool],1,2);
            d(d<0) = 0;
            stat1 = sum(d,2)./f; 

            pool = find(~ismember(1:n,S1PCs));
        case 'si'
            load('SpatialInfo.mat','MI');
            n = length(MI);
            stat1 = MI; 
            
            pool = find(~ismember(1:n,S1PCs));
        case 'ti'
            load('TemporalInfo.mat','MI');
            n = length(MI);
            stat1 = MI; 
      
            pool = find(~ismember(1:n,S1TCs));
    end
    
    %Z-score the statistic relative to the pool of neurons in session 1
    %that are currently not place cells. This will include the cells that
    %are not categorized as place cells but eventually become categorizd as
    %place cells in the subsequent session.
    %stat(pool) = zscore(stat(pool));
    newPCStatS1 = stat1(newPCsS1);   
    
    %Randomly sample from the pool of all cells. 
    good = EliminateUncertainMatches([base,comp],1:n);
    r1 = randsample(stat1(good),1000,true);
    
%     load('TemporalInfo.mat','MI');
%     MI(pool) = zscore(MI(pool)); 
%     r = MI(newPCs);
    
    [~,p] = kstest2(r1,newPCStatS1);
    
    if plotit
        figure; hold on
        histogram(r1,'normalization','probability','binwidth',binsize,'edgecolor','none');
        histogram(newPCStatS1,'normalization','probability','binwidth',binsize,'edgecolor','none');
        title(['P = ',num2str(p)]);
    end
    
%% Get stats of new time cells. 
    cd(comp.Location);
    switch statType
        case 'fr'
            load('Pos_align.mat','PSAbool');
            [n,f] = size(PSAbool);
            d = diff([zeros(n,1) PSAbool],1,2);
            d(d<0) = 0;
            stat2 = sum(d,2)./f; 

            pool = find(~ismember(1:n,S2PCs));
        case 'ti'
            load('TemporalInfo.mat','MI');
            stat2 = MI; 
            n = length(MI);
            
            pool = find(~ismember(1:n,S2TCs));            
        case 'si'
            load('SpatialInfo.mat','MI');
            stat2 = MI; 
            
            n = length(MI);
            pool = find(~ismember(1:n,S2PCs)); 
    end
    
    good = EliminateUncertainMatches([comp,base],1:n);
    r2 = randsample(stat2(good),1000,true);
    newPCStatS2 = stat2(newPCsS2);
    
    
end