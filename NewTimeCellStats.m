function [newTCStatS1,r1,stat1,newTCStatS2,r2,stat2] = NewTimeCellStats(base,comp,statType,varargin)
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
    [newTCsS1,newTCsS2] = getNewTimeCells(base,comp);
    S1TCs = getTimeCells(base);
    S2TCs = getTimeCells(comp);
    S1PCs = getPlaceCells(base,.01);
    S2PCs = getPlaceCells(comp,.01);
    cd(base.Location);

    binsize = .1;
    switch statType
        case 'fr'
            load('Pos_align.mat','PSAbool');
            [n,f] = size(PSAbool);
            d = diff([zeros(n,1) PSAbool],1,2);
            d(d<0) = 0;
            stat1 = sum(d,2)./f; 

            pool = find(~ismember(1:n,S1TCs));
        case 'ti'
            load('TemporalInfo.mat','MI');
            n = length(MI);
            stat1 = MI; 
            
            pool = find(~ismember(1:n,S1TCs));
        case 'si'
            load('SpatialInfo.mat','MI');
            n = length(MI);
            stat1 = MI; 

            pool = find(~ismember(1:n,S1PCs));
    end
    
    %Z-score the statistic relative to the pool of neurons in session 1
    %that are currently not place cells. This will include the cells that
    %are not categorized as place cells but eventually become categorizd as
    %place cells in the subsequent session.
    %stat(pool) = zscore(stat(pool));
    newTCStatS1 = stat1(newTCsS1);   
    
    %Randomly sample from the pool of all cells. 
    good = EliminateUncertainMatches([base,comp],1:n);
    %good = intersect(good,S1TCs);
    r1 = randsample(stat1(good),1000,true);
    
%     load('SpatialInfo.mat','MI');
%     MI(pool) = zscore(MI(pool)); 
%     r = MI(newTCs);

    if plotit
        [~,p] = kstest2(r1,newTCStatS1);
        figure; hold on
        histogram(r1,'normalization','probability','binwidth',binsize,'edgecolor','none');
        histogram(newTCStatS1,'normalization','probability','binwidth',binsize,'edgecolor','none');
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

            pool = find(~ismember(1:n,S2TCs));
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
    %good = intersect(good,S2TCs);
    r2 = randsample(stat2(good),1000,true);
    newTCStatS2 = stat2(newTCsS2);
    
    
end