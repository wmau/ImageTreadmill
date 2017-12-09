function [r,pval,p] = corrInfo(mds,cellType)
%
%
%

%%
    DATA = CompileMultiSessionData(mds,{'si','ti','timecells','placecells'}); 
    B = 1000;
    corrtype = 'pearson';
    
    switch cellType
        case {'dual','either'}, c = 'k';
        case 'time', c = [0 .5 .5];
        case 'place', c = [.58 .44 .86]; 
    end
    
    nSessions = length(mds);
    neurons = cell(1,nSessions);
    figure; hold on;
    for s = 1:nSessions
        switch cellType
            case 'dual', neurons{s} = intersect(DATA.timecells{s},DATA.placecells{s});
            case 'time', neurons{s} = setdiff(DATA.timecells{s},DATA.placecells{s});
            case 'place', neurons{s} = setdiff(DATA.placecells{s},DATA.timecells{s});
            case 'either', neurons{s} = union(DATA.timecells{s},DATA.placecells{s});
        end
    end
    
    for s = 1:nSessions
        DATA.si{s}(neurons{s}) = zscore(DATA.si{s}(neurons{s}));
        DATA.ti{s}(neurons{s}) = zscore(DATA.ti{s}(neurons{s}));
    end
    
    SI = cell2mat(cellfun(@(x,y) x(y),DATA.si,neurons,'unif',0)');
    TI = cell2mat(cellfun(@(x,y) x(y),DATA.ti,neurons,'unif',0)');
    
    h = scatter(SI,TI,50,c);
    set(gca,'tickdir','out','fontsize',12,'linewidth',4);
    alpha(h,.5); 
    lsline;
    xlabel('Norm. Spatial Info','fontsize',15);
    ylabel('Norm. Temporal Info','fontsize',15);
    
    [r,pval] = corr(SI,TI,'type',corrtype);
    title(['R = ',num2str(r),', p = ',num2str(pval)],'fontsize',15);
%     r = r^2
%     pval 
    
    shuffle = zeros(1,B); 
    for i=1:B
        shuffle(i) = corr(SI,TI(randperm(length(TI))),'type',corrtype);
    end
%     shuffle = shuffle.^2;
    p = sum(r<shuffle)/B;
    
end