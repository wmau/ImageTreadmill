function [newTCStat,r] = NewTimeCellStats(base,comp,statType,varargin)
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
    newTCs = getNewTimeCells(base,comp);
    S1TCs = getTimeCells(base);
    S1PCs = getPlaceCells(base,.01);

    binsize = .1;
    switch statType
        case 'fr'
            load('Pos_align.mat','PSAbool');
            [n,f] = size(PSAbool);
            d = diff([zeros(n,1) PSAbool],1,2);
            d(d<0) = 0;
            stat = sum(d,2)./f; 

            pool = find(~ismember(1:n,S1TCs));
        case 'ti'
            load('TemporalInfo.mat','MI');
            n = length(MI);
            stat = MI; 
            
            pool = find(~ismember(1:n,S1TCs));
        case 'si'
            load('SpatialInfo.mat','MI');
            n = length(MI);
            stat = MI; 

            pool = find(~ismember(1:n,S1PCs));
    end
    
    stat(pool) = zscore(stat(pool));
    r = randsample(stat(pool),1000,true);
    newTCStat = stat(newTCs);   
    [~,p] = kstest2(r,newTCStat);
    
    if plotit
        figure; hold on
        histogram(r,'normalization','probability','binwidth',binsize,'edgecolor','none');
        histogram(newTCStat,'normalization','probability','binwidth',binsize,'edgecolor','none');
        title(['P = ',num2str(p)]);
    end
end