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
    
%% Determine what type of neuron to analyze.
    switch stabilityType
        case 'time'
            switch statType
                case 'ti',cellGet = 'timecells'; 
                case 'si',cellGet = 'timecells'; 
                case {'fr','fluor'}, cellGet = 'timecells';
            end
        case 'place'
            switch statType
                case 'ti',cellGet = 'placecells'; 
                case 'si',cellGet = 'placecells'; 
                case {'fr','fluor'}, cellGet = 'placecells';
            end
    end
   
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
 
            neurons{a}{s} = AcquireTimePlaceCells(mds(ssns(s)),cellGet);   
            neurons{a}{s} = EliminateUncertainMatches([mds(ssns(s)) mds(ssns(s+1))],neurons{a}{s});
                             
            %% Get information scores and neurons. 
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI');
                    stats{a}{s} = MI; 
                    
                case 'si'
                    load('SpatialInfo.mat','MI');
                    stats{a}{s} = MI;       
                    
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [n,f] = size(PSAbool);
                    d = diff([zeros(n,1) PSAbool],1,2);
                    d(d<0) = 0;
                    stats{a}{s} = sum(d,2)./f;
                    
            end  
            
            %% Correlate tuning curves. 
            switch stabilityType
                case 'time'
                    corrs{a}{s} = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons{a}{s});
                    
                case 'place'
                    corrs{a}{s} = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons{a}{s});
            end
            
            %% Normalize.
%             if strcmp(statType,'fr')
%                 stats{a}{s} = (stats{a}{s}-min(stats{a}{s}))./range(stats{a}{s}); 
%             else
                stats{a}{s}(neurons{a}{s}) = zscore(stats{a}{s}(neurons{a}{s}));
%            end

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