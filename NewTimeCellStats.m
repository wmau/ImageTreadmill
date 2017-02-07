function [newTCStat,r] = NewTimeCellStats(base,comp,statType)
%
%
%

%%
    cd(base.Location);
    newTCs = getNewTimeCells(base,comp);
    oldTCs = getTimeCells(base);
    
    [~,TMactive] = nNeuronsActiveonTM(base);
    TMactive = setdiff(TMactive,oldTCs);
    
    switch statType
        case 'fr'
            load('Pos_align.mat','PSAbool');
            [n,f] = size(PSAbool);
            d = diff([zeros(n,1) PSAbool],1,2);
            d(d<0) = 0;
            stat = sum(d,2)./f; 
        case 'ti'
            load('TemporalInfo.mat','MI');
            stat = MI; 
    end
    
    newTCsandActiveTM = union(TMactive,newTCs);
    stat(newTCsandActiveTM) = zscore(stat(newTCsandActiveTM));
    
    r = randsample(stat(TMactive),1000,true);
    newTCStat = stat(newTCs);   
    
end