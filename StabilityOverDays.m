function StabilityOverDays(mds)
%
%
%

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals); 
    B = 100;
    
    [corrs,shuffleCorrs] = deal(cell(nAnimals,1));
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions = length(ssns);
        
        [corrs{a},shuffleCorrs{a}] = deal(cell(nSessions-1,1));
        TCs = getTimeCells(mds(ssns(1)));
        
        for s=2:nSessions
            disp(['Analyzing ',mds(ssns(s)).Animal,' on ',mds(ssns(s)).Date,'.']);
            cd(mds(ssns(s)).Location);
                       
            TCs = EliminateUncertainMatches([mds(ssns(1)) mds(ssns(s))],TCs);
            
            temp = CorrTrdmllTrace(mds(ssns(1)),mds(ssns(s)),TCs);
            corrs{a}{s-1} = temp(:,1);
            
            %Preallocate.
            shuffleCorrs{a}{s-1} = nan(size(temp,1),B);
            
            %Progress bar. 
            resolution = 2;
            updateInc = round(B/(100/resolution));
            p = ProgressBar(100/resolution);
            for i=1:B
                temp = CorrTrdmllTrace(mds(ssns(1)),mds(ssns(s)),TCs,'shuffle',true);
                shuffleCorrs{a}{s-1}(:,i) = temp(:,1);
                
                if round(i/updateInc) == (i/updateInc), p.progress; end
            end
        end
        
        %Mean of all surrogate correlations. Rows = sessions,
        meanShuffleCorrs{a} = cell2mat(cellfun(@(x) nanmean(x),shuffleCorrs{a},'unif',0));
        meanCorrs{a} = cellfun(@(x) x(1), cellfun(@(x) nanmean(x),corrs{a},'unif',0));
    end
    
%% 
    %Reorganize the mean shuffled correlations. 
    daysOut = max(cellfun(@(x) x(1), cellfun(@size, meanShuffleCorrs,'unif',0)));
    [r,m,color,p] = deal(cell(daysOut,1));
    c = parula(nAnimals); 
    for a=1:nAnimals
        for d=1:size(meanShuffleCorrs{a},1)
            r{d} = [r{d}, meanShuffleCorrs{a}(d,:)];
            m{d} = [m{d}, meanCorrs{a}(d)];
            color{d} = [color{d}; c(a,:)];
            %p{d} = [p{d}, sum(meanCorrs{a}(d)<meanShuffleCorrs{a}(d,:))/size(meanShuffleCorrs{a},2)];
        end
    end
            
    %Get mean of all correlation means.
    M = cellfun(@mean,m);
    
    %Sort. 
    sortedR = cellfun(@sort,r,'unif',0); 
    Mr = cellfun(@mean,r);
    %Lerr = Mr - cellfun(@(x) x(round(.05*length(x))),sortedR);
    %Uerr = cellfun(@(x) x(round(.95*length(x))),sortedR) - Mr;
    
    figure('Position',[610 270 400 350]); hold on;
    b = bar([1:daysOut],M,'facecolor','none','linewidth',2,'edgecolor',[0 .5 .5]);
    w = b.BarWidth;
    for d=1:daysOut
        line([d-w/2 d+w/2],[Mr(d) Mr(d)],'color','b','linewidth',2);
        errorbar(d,M(d),std(m{d})/sqrt(length(m{d})),'linewidth',2,'color',[0 .5 .5]);
        scat = scatter(d*ones(length(m{d}),1),m{d}',[],color{d},'d','filled');
        alpha(scat,.5);
    end
    set(gca,'xtick',[1:daysOut],'tickdir','out');
    xlabel('Days from Reference'); 
    ylabel('Mean Corr. Coeff.');
    set(gca,'linewidth',2)
end