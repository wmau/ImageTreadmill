function [STATS,nNeurons,stable,unstable] = PartitionStats(mds,stabilityType,statType)
%[STATS,nNeurons,stable,unstable] = PartitionStats(mds,celltype,infotype)
%
%   Partitions an arbitrary neuron statistic based on stability in the
%   domain specified by stabilityCriterion.
%
%   INPUTS
%       mds: Session entries which you want the information content for.
%
%       celltype: 'time' or 'space' to select which criterion you want to
%       use for stability. 
%
%       infotype: Information in either the 'time' or 'space' dimension.
%
%   OUTPUTS
%       STATS: Struct with fields stable or unstable each being a cell array
%       with entries for each session containing the information content of
%       cells that are either stable or unstable. 
%
%       nNeurons:
%
%       stable & unstable: 

%% Set up.
    animals = unique({mds.Animal});
    nAnimals = length(animals); 
    statType = lower(statType);

%% Compile
    [STATS.stable,STATS.unstable,stable,unstable] = deal(cell(1,nAnimals));
    [nNeurons.stable,nNeurons.unstable] = deal(zeros(1,nAnimals)); 
    
    %p-value criterion to be considered a place cell. 
    PCcrit = .01;
    for a = 1:nAnimals       
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
    
        [nStable,nUnstable] = deal(0); 
        [STATS.stable{a},STATS.unstable{a}] = deal([]);
        [stable{a},unstable{a}] = deal(cell(1,length(ssns)-1));
        for s = 1:length(ssns)-1
            cd(mds(ssns(s)).Location);
 
            %Get place and time cells.
            PCs = getPlaceCells(mds(ssns(s)),PCcrit);
            TCs = getTimeCells(mds(ssns(s)));
                        
            switch statType
                %Temporal information. 
                case 'ti' 
                    load('TemporalInfo.mat','MI','Ispk','Isec');
                    stat = MI;                

                    %If looking at temporal information and time stability,
                    %get all time cells. Otherwise, get dual time/place
                    %cells. 
                    switch stabilityType 
                        case 'time',neurons = TCs;
                        case 'place',neurons = intersect(TCs,PCs);
                    end
                    
                %Spatial information.
                case 'si'
                    load('SpatialInfo.mat','MI','Ispk','Isec');
                    stat = MI;

                    %If looking at spatial information and spatial
                    %stability, get all place cells. Otherwise, get dual
                    %time/place cells. 
                    switch stabilityType 
                        case 'time',neurons = intersect(TCs,PCs);
                        case 'place',neurons = PCs;
                    end
                
                %Firing rate (really Ca event rate).
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [n,f] = size(PSAbool);
                    d = diff([zeros(n,1) PSAbool],1,2);
                    d(d<0) = 0;
                    stat = sum(d,2)./f; 
                    
                    %Select cells based on stability type. 
                    switch stabilityType
                        case 'time', neurons = TCs;
                        case 'place',neurons = PCs;
                    end
                    
                %Mean fluorescence.
                case 'fluor'
                    load('Pos_align.mat','DFDTtrace');
                    stat = mean(DFDTtrace,2);
                
                    %Select cells based on stability type. 
                    switch stabilityType
                        case 'time', neurons = TCs;
                        case 'place',neurons = PCs;
                    end
                    
                %Time peak. 
                case 'pk'
                    [~,stat] = getTimePeak(mds(ssns(s)));

                    %Should only ever be time cells. 
                    neurons = TCs;                              
            end
                     
            neurons = EliminateUncertainMatches([mds(ssns(s)) mds(ssns(s+1))],neurons);
            
            %stat = (stat-min(stat))./range(stat);
            %Normalize.
            if strcmp(statType,'fr')
                %stat = (stat-min(stat))./range(stat);
                stat(neurons) = (stat(neurons)-min(stat(neurons)))./range(stat(neurons));
            else
                stat(neurons) = zscore(stat(neurons));
            end
                           
            switch stabilityType
                case 'time'           
                    %Get correlation coefficients and p-values. 
                    corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons);

                case 'place'
                    %Get the correlation coefficients and p-values.
                    corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons);       
            end
            
            %[~,stblcrit] = fdr_bh(corrStats(~isnan(corrStats(:,2)),2));
            stblcrit = .01/length(neurons);

            %Stable time cells based on correlation.
            stable{a}{s} = intersect(find(corrs(:,2) < stblcrit),neurons);
            unstable{a}{s} = intersect(find(corrs(:,2) >= stblcrit | isnan(corrs(:,2))),neurons);
                 
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
    end
    
end
        
        