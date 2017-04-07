function StabilityOverDays(mds)
%StabilityOverDays(mds)
%
%   Plots the mean correlation coefficient of pairwise time cell ensemble
%   correlations. 

%% 
    animals = unique({mds.Animal});     %Cell array of animal names. 
    nAnimals = length(animals); 
    B = 100;
    teal = [0 .5 .5];
    
    [corrs,pval,stability,shuffleCorrs,meanShuffleCorrs,meanCorrs] = deal(cell(nAnimals,1));
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions = length(ssns);
        
        %Preallocate. 
        [corrs{a},pval{a},stability{a},shuffleCorrs{a}] = deal(cell(nSessions-1,1));
        TCs = getTimeCells(mds(ssns(1)));
        %TCs = getPlaceCells(mds(ssns(1)),.01);
        
        for s=2:nSessions
            disp(['Analyzing ',mds(ssns(s)).Animal,' on ',mds(ssns(s)).Date,'.']);
            cd(mds(ssns(s)).Location);
                           
            %Get rid of unmapped neurons. 
            TCs = EliminateUncertainMatches([mds(ssns(1)) mds(ssns(s))],TCs);
            
            %Critical p-value to be considered stable. 
            stableCrit = 0.01/length(TCs);
            
            %Correlation coefficients and p-values.
            temp = CorrTrdmllTrace(mds(ssns(1)),mds(ssns(s)),TCs);
            corrs{a}{s-1} = temp(:,1);
            pval{a}{s-1} = temp(:,2); 
            stability{a}{s-1} = temp(:,2)<stableCrit;
            
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
    [r,m,corrDays,stabilityDays,color,p] = deal(cell(daysOut,1));
    c = parula(nAnimals); 
    for a=1:nAnimals
        for d=1:size(meanShuffleCorrs{a},1)
            r{d} = [r{d}, meanShuffleCorrs{a}(d,:)];    %Mean corr coeff for each shuffle iteration.
            m{d} = [m{d}, meanCorrs{a}(d)];             %Real mean corr coeff.
            corrDays{d} = [corrDays{d}; corrs{a}{d}];    %All corr coeffs. 
            %Stability labels. 
            stabilityDays{d} = [stabilityDays{d}; stability{a}{d}];
            color{d} = [color{d}; c(a,:)];              %Color for labeling animals. 
            %p-value each day. 
            p{d} = [p{d}, sum(meanCorrs{a}(d)<meanShuffleCorrs{a}(d,:))/size(meanShuffleCorrs{a},2)];
        end
    end
            
    %Turn into logicals. 
    stabilityDays = cellfun(@logical,stabilityDays,'unif',false);
    %Get mean of all correlation means.
    M = cellfun(@mean,m);
    
    %Sort. 
    sortedR = cellfun(@sort,r,'unif',0); 
    Mr = cellfun(@mean,r);
    %Lerr = Mr - cellfun(@(x) x(round(.05*length(x))),sortedR);
    %Uerr = cellfun(@(x) x(round(.95*length(x))),sortedR) - Mr;
    
    figure('Position',[610 270 400 350]); hold on;
    b = bar([1:daysOut],M,'facecolor','none','linewidth',2,'edgecolor',teal);
    w = b.BarWidth;
    for d=1:daysOut
        line([d-w/2 d+w/2],[Mr(d) Mr(d)],'color','b','linewidth',2);
        errorbar(d,M(d),std(m{d})/sqrt(length(m{d})),'linewidth',2,'color',teal);
        
        jitter = d - (.1*randn(length(corrDays{d}),1));
        s(1) = scatter(jitter(stabilityDays{d}),corrDays{d}(stabilityDays{d}),20,teal,'filled');
        s(2) = scatter(jitter(~stabilityDays{d}),corrDays{d}(~stabilityDays{d}),20,'r','x');
        
        s(3) = scatter(d*ones(length(m{d}),1),m{d}',120,color{d},'d','filled');
        alpha(s,.5);
    end
    b = bar([1:daysOut],M,'facecolor','none','linewidth',2,'edgecolor',teal);
    set(gca,'xtick',[1:daysOut],'tickdir','out');
    xlabel('Days from Reference'); 
    ylabel('Mean Corr. Coeff.');
    set(gca,'linewidth',4,'fontsize',15);
end