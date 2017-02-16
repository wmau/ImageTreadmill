function [STATS,nNeurons,stable,unstable] = StabilityStatsEarlyLate(mds,stabilityType,statType,treadmillTime)
%
%
%

%% Set up.
    animals = unique({mds.Animal});
    nAnimals = length(animals);
    stabilityType = lower(stabilityType); 
    statType = lower(statType);

%% Determine what type of neuron to analyze.
    switch statType
        case 'ti',cellGet = 'timecells'; 
        case 'si',cellGet = 'timecells'; 
        case {'fr','fluor'}, cellGet = 'timecells';
    end
   
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
       
            %Get time and place cells. 
            TCs = getTimeCells(mds(ssns(s)));
 
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI');
                    stat = MI; 
                                  
                case 'si'
                    load('SpatialInfo.mat','MI');
                    stat = MI;

                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [n,f] = size(PSAbool);
                    d = diff([zeros(n,1) PSAbool],1,2);
                    d(d<0) = 0;
                    stat = sum(d,2)./f;  
            end
            
            neurons = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
            neurons = EliminateUncertainMatches([mds(ssns(s)) mds(ssns(s+1))],neurons);
            
            %Get early/late.
            [~,order] = PastalkovaPlot(mds(ssns(s)),'plotit',false);
            middle = length(order)/2;

            if strcmp(treadmillTime,'early')
                neurons = intersect(neurons,TCs(order < middle));
            elseif strcmp(treadmillTime,'late')
                neurons = intersect(neurons,TCs(order >= middle)); 
            end
 
            %Normalize.
%             if strcmp(statType,'fr')
%                 stat = (stat-min(stat))./range(stat);
%             else
                stat(TCs) = zscore(stat(TCs));
%             end

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