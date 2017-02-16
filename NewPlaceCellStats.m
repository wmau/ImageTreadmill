function [newPCStat,r] = NewPlaceCellStats(base,comp,statType,varargin)
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
    cd(base.Location);
    newPCs = getNewPlaceCells(base,comp);
    S1PCs = getPlaceCells(base,.01);
    S1TCs = getTimeCells(base);

    binsize = .1;
    switch statType
        case 'fr'
            load('Pos_align.mat','PSAbool');
            [n,f] = size(PSAbool);
            d = diff([zeros(n,1) PSAbool],1,2);
            d(d<0) = 0;
            stat = sum(d,2)./f; 

            pool = find(~ismember(1:n,S1PCs));
        case 'si'
            load('SpatialInfo.mat','MI');
            n = length(MI);
            stat = MI; 
            
            pool = find(~ismember(1:n,S1PCs));
        case 'ti'
            load('TemporalInfo.mat','MI');
            n = length(MI);
            stat = MI; 

            %If looking cross modally, look at all neurons.
            pool = find(~ismember(1:n,S1TCs));
    end
    
    %Z-score the statistic relative to the pool of neurons in session 1
    %that are currently not place cells. This will include the cells that
    %are not categorized as place cells but eventually become categorizd as
    %place cells in the subsequent session.
    stat(pool) = zscore(stat(pool));
    
    %Randomly sample from the pool of not-place cells. 
    r = randsample(stat(pool),1000,true);
    newPCStat = stat(newPCs);   
    [~,p] = kstest2(r,newPCStat);
    
    if plotit
        figure; hold on
        histogram(r,'normalization','probability','binwidth',binsize,'edgecolor','none');
        histogram(newPCStat,'normalization','probability','binwidth',binsize,'edgecolor','none');
        title(['P = ',num2str(p)]);
    end
end