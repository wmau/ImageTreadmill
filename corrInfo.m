function corrInfo(mds)
%
%
%

%%
    DATA = CompileMultiSessionData(mds,{'si','ti','timecells','placecells'}); 
    B = 10000;
    corrtype = 'spearman';
    
    nSessions = length(mds);
    both = cell(1,nSessions);
    PCs = cell(1,nSessions);
    TCs = cell(1,nSessions); 
    either = cell(1,nSessions);
    figure; hold on;
    for s = 1:nSessions
        both{s} = intersect(DATA.timecells{s},DATA.placecells{s});
        PCs{s} = setdiff(DATA.placecells{s},DATA.timecells{s});
        TCs{s} = setdiff(DATA.timecells{s},DATA.placecells{s});
        either{s} = union(DATA.timecells{s},DATA.placecells{s});   
    end
    
    SI = cell2mat(cellfun(@(x,y) x(y),DATA.si,both,'unif',0)');
    TI = cell2mat(cellfun(@(x,y) x(y),DATA.ti,both,'unif',0)');
    
    h = scatter(SI,TI,10,'k','filled');
    alpha(h,.5); 
    
    [r,pval] = corr(SI,TI,'type',corrtype); 
    r = r^2
    pval 
    
    shuffle = zeros(1,B); 
    for i=1:B
        shuffle(i) = corr(SI,TI(randperm(length(TI))),'type',corrtype);
    end
    shuffle = shuffle.^2;
    p = sum(r<shuffle)/B
    
end