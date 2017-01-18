function CORRS = corrInfoStability(mds,stabilityType,statType)
%corrInfoStability(mds,stabilityType,statType)
%
%   Correlate the correlation coefficient of a cell and its temporal or
%   spatial information. Takes same inputs as PartitionStats. Work on
%   commenting this.

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals);
    statType = lower(statType);
    stabilityType = lower(stabilityType); 
    PCcrit = .01;
    
%% 
    [stats,corrs,neurons] = deal(cell(nAnimals,1)); 
    STATS = []; 
    CORRS = [];
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal}));
        
        [stats{a},corrs{a},neurons{a}] = deal(cell(length(ssns)-1,1));
        for s=1:length(ssns)-1
            cd(mds(ssns(s)).Location);
            
            load('TimeCells.mat','TimeCells');
            load('TemporalInfo.mat','sig'); 
            load('Placefields.mat','pval');
            load('PlacefieldStats.mat','PFnHits','PFpcthits','bestPF');
            load('SpatialInfo.mat','MI');
            
            %Matrix index. 
            idx = sub2ind(size(PFnHits),1:size(PFnHits,1),bestPF');
            
            %Get time and place cells. 
            TCs = intersect(TimeCells,find(sig)); 
            PCs = find(pval<PCcrit & MI'>0 & PFnHits(idx)>4); 
            
            %% Get information scores and neurons. 
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI');
                    stats{a}{s} = MI; 
                    
                    switch stabilityType
                        case 'time'
                            neurons{a}{s} = TCs;
                        case 'place'
                            neurons{a}{s} = intersect(TCs,PCs);
                    end
                    
                case 'si'
                    load('SpatialInfo.mat','MI');
                    stats{a}{s} = MI;
                    
                    switch stabilityType 
                        case 'time'
                            neurons{a}{s} = intersect(TCs,PCs); 
                        case 'place'
                            neurons{a}{s} = PCs'; 
                    end
                    
                case 'fr'
                    load('Pos_align.mat','PSAbool')
                    [~,f] = size(PSAbool);
                    
                    stats{a}{s} = sum(PSAbool,2)./f;
                    
                    switch stabilityType
                        case 'time'
                            neurons{a}{s} = TCs;
                        case 'place'
                            neurons{a}{s} = PCs'; 
                    end
            end
            
            %% Correlate tuning curves. 
            switch stabilityType
                case 'time'
                    corrs{a}{s} = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons{a}{s});
                    
                case 'place'
                    corrs{a}{s} = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons{a}{s});
            end
            
            %% Normalize.
            if strcmp(statType,'fr')
                stats{a}{s} = (stats{a}{s}-min(stats{a}{s}))./range(stats{a}{s}); 
            else
                stats{a}{s}(neurons{a}{s}) = zscore(stats{a}{s}(neurons{a}{s}));
            end

        end
        
        STATS = [STATS; cell2mat(cellfun(@(x,y) x(y),stats{a},neurons{a},'unif',0))];
        CORRS = [CORRS; cell2mat(cellfun(@(x,y) x(y),corrs{a},neurons{a},'unif',0))];
    end
    
    CORRS(isnan(CORRS)) = 0;
    
    [pr,pp] = corr(STATS,CORRS,'rows','complete','type','pearson');
    [sr,sp] = corr(STATS,CORRS,'rows','complete','type','spearman');
    [kr,kp] = corr(STATS,CORRS,'rows','complete','type','kendall');

    figure('Name',[statType, ' vs ',stabilityType, ' stability']);
    scatter(STATS,CORRS,'.');
    ylabel('Correlation Rho');
    xlabel('Information [z-scored bits]');
    title({['Pearson R = ',num2str(pr),', p = ',num2str(pp)],...
           ['Spearman R = ',num2str(sr),', p = ',num2str(sp)],...
           ['Kendall R = ',num2str(kr),', p = ',num2str(kp)]});
end