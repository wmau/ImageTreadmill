function [STATS,nNeurons,stable,unstable] = StabilityStatsEarlyLate(mds,stabilityType,statType,treadmillTime)
%
%
%

%% Set up.
    animals = unique({mds.Animal});
    nAnimals = length(animals);
    stabilityType = lower(stabilityType); 
    statType = lower(statType);

%% 
    [STATS.stable,STATS.unstable,stable,unstable] = deal(cell(1,nAnimals));
    [nNeurons.stable,nNeurons.unstable] = deal(zeros(1,nAnimals)); 
    
    for a=1:nAnimals
        %Get all sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal}));
        nSessions = length(ssns);
        
        [nStable,nUnstable] = deal(0);
        [STATS.stable{a},STATS.unstable{a}] = deal([]);
        [stable{a},unstable{a}] = deal(cell(1,nSessions-1));
        
        for s=1:nSessions-1
            cd(mds(ssns(s)).Location);
                     
            %Critical p-value for place cell consideration.
            PCcrit = .01;
            
            %Get time and place cells. 
            TCs = getTimeCells(mds(ssns(s)));
            PCs = getPlaceCells(mds(ssns(s)),PCcrit);
            
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI');
                    stat = MI; 
                    
                    if strcmp(stabilityType,'time')
                        neurons = TCs; 
                    elseif strcmp(stabilityType,'place')
                        neurons = intersect(TCs,PCs);
                    end
                    
                case 'si'
                    load('SpatialInfo.mat','MI');
                    stat = MI;
                    
                    if strcmp(stabilityType,'place')
                        neurons = PCs;
                    elseif strcmp(stabilityType,'time')
                        neurons = intersect(TCs,PCs);
                    end
                    
                case 'pk'
                    [~,stat] = getTimePeak(mds(ssns(s))); 
                    
                    neurons = TCs;
                    
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [~,f] = size(PSAbool); 
                    stat = sum(PSAbool,2)./f;
                    
                    if strcmp(stabilityType,'time')
                        neurons = TCs;
                    elseif strcmp(stabilityType,'place')
                        neurons = PCs;
                    end
            end
            
            if strcmp(stabilityType,'time') || ...
              (strcmp(stabilityType,'place') && strcmp(statType,'time'))
                %Get early/late.
                [~,order] = PastalkovaPlot(mds(ssns(s)),'plotit',false);
                middle = length(order)/2;

                if strcmp(treadmillTime,'early')
                    neurons = intersect(neurons,TCs(order < middle));
                elseif strcmp(treadmillTime,'late')
                    neurons = intersect(neurons,TCs(order >= middle)); 
                end
            end
            
            %Normalize.
            if strcmp(statType,'fr')
                stat = (stat-min(stat))./range(stat);
            else
                stat(neurons) = zscore(stat(neurons));
            end
            
            
            if strcmp(stabilityType,'time')           
                %Get correlation coefficients and p-values. 
                corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons);
                
                %[~,stblcrit] = fdr_bh(corrStats(~isnan(corrStats(:,2)),2));
                stblcrit = .01/length(neurons);

                %Stable time cells based on correlation and non-shifting time
                %field.
                stable{a}{s} = intersect(find(corrs(:,2) < stblcrit),neurons);
                unstable{a}{s} = intersect(find(corrs(:,2) >= stblcrit | isnan(corrs(:,2))),neurons);
         
            elseif strcmp(stabilityType,'place')
                
                %Get the correlation coefficients and p-values.
                corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons);
                
                %[~,stblcrit] = fdr_bh(corrStats(~isnan(corrStats(:,2)),2));
                stblcrit = .01/length(neurons);
                
                %Stable place cells based on correlation p-value.
                stable{a}{s} = intersect(find(corrs(:,2) < stblcrit),neurons); 
                unstable{a}{s} = intersect(find(corrs(:,2) >= stblcrit | isnan(corrs(:,2))),neurons);     
            end
           
            %Get number of stable and unstable cells. Add this one per
            %session but track the number per animal. 
            nStable = nStable + length(stable{a}{s}); 
            nUnstable = nUnstable + length(unstable{a}{s}); 
            
            %Get the temporal information values. 
            STATS.stable{a} = [STATS.stable{a}; stat(stable{a}{s})];
            STATS.unstable{a} = [STATS.unstable{a}; stat(unstable{a}{s})];
        end
        
        %Get number of neurons that are (un)stable per animal. 
        nNeurons.stable(a) = nStable; 
        nNeurons.unstable(a) = nUnstable;
        
        %After one animal.
    end
            
end