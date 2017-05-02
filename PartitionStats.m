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
    
%% Determine what type of neuron to analyze.
    switch stabilityType
        case 'time'
            switch statType
                case 'ti',cellGet = 'timecells'; 
                case 'si',cellGet = 'timecells'; 
                case {'fr','fluor','tfw','pfw'}, cellGet = 'timecells';
            end
        case 'place'
            switch statType
                case 'ti',cellGet = 'placecells'; 
                case 'si',cellGet = 'placecells'; 
                case {'fr','fluor','tfw','pfw'}, cellGet = 'placecells';
            end
    end
%% Compile
    [STATS.stable,STATS.unstable,stable,unstable] = deal(cell(1,nAnimals));
    [nNeurons.stable,nNeurons.unstable] = deal(zeros(1,nAnimals)); 
    
    %p-value criterion to be considered a place cell. 
    for a = 1:nAnimals       
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
    
        [nStable,nUnstable] = deal(0); 
        [STATS.stable{a},STATS.unstable{a}] = deal([]);
        [stable{a},unstable{a}] = deal(cell(1,length(ssns)-1));
        for s = 1:length(ssns)-1
            cd(mds(ssns(s)).Location);
      
            neurons = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
            neurons = EliminateUncertainMatches([mds(ssns(s)) mds(ssns(s+1))],neurons);  
                        
            switch statType
                %Temporal information. 
                case 'ti' 
                    load('TemporalInfo.mat','MI','Ispk','Isec');
                    stat = MI; 
            
                %Spatial information.
                case 'si'
                    load('SpatialInfo.mat','MI','Ispk','Isec');
                    stat = MI;
                    
                %Time field width.
                case 'tfw'
%                     load('TimeCells','curves');
%                     stat = cellfun(@sum, cellfun(@(x) x > 0, curves.tuning,'unif',0));
                    load('TimefieldStats.mat','TFpcthits');
                    stat = TFpcthits;

                %Place field width. 
                case 'pfw' 
                    load('Placefields.mat','TMap_gauss');
                    load('PlacefieldStats.mat','PFnHits','PFarea','PFpcthits');
                    [~,bestPF] = max(PFnHits,[],2);
                    idx = sub2ind(size(PFnHits),1:size(PFnHits,1),bestPF');
                    %stat = cell2mat(cellfun(@(x) sum(x(:) > 0), TMap_gauss,'unif',0))';
                    stat = PFpcthits(idx)';
                    stat(isnan(stat)) = 0;
                    
                %Firing rate (really Ca event rate).
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [n,f] = size(PSAbool);
                    d = diff([zeros(n,1) PSAbool],1,2);
                    d(d<0) = 0;
                    stat = sum(d,2)./f;                  
                    
                %Mean fluorescence.
                case 'fluor'
                    load('Pos_align.mat','DFDTtrace');
                    stat = mean(DFDTtrace,2);
                                
                %Time peak. 
                case 'pk'
                    [~,stat] = getTimePeak(mds(ssns(s)));
                    
                    %Get time cells.
                    TCs = getTimeCells(mds(ssns(s)));

                    %Should only ever be time cells. 
                    neurons = TCs;                              
            end
                     
            
            %Normalize.
            %if strcmp(statType,'fr')
                %stat = (stat-min(stat))./range(stat);
                %stat(neurons) = (stat(neurons)-min(stat(neurons)))./range(stat(neurons));
            %else
%             if any(strcmp(statType,{'tfw','pfw'}))
%                 switch stabilityType
%                     case 'time'
%                         switch statType
%                             case 'tfw', stat(neurons) = zscore(stat(neurons));
%                             case 'pfw'
%                                 PCs = getPlaceCells(mds(ssns(s)),.01);
%                                 mu = mean(stat(PCs));
%                                 sigma = std(stat(PCs));
%                                 stat(neurons) = bsxfun(@rdivide, bsxfun(@minus, stat(neurons), mu),sigma);
%                         end
%                     case 'place'
%                         switch statType
%                             case 'tfw'
%                                 TCs = getTimeCells(mds(ssns(s)));
%                                 mu = mean(stat(TCs));
%                                 sigma = std(stat(TCs));
%                                 stat(neurons) = bsxfun(@rdivide, bsxfun(@minus, stat(neurons), mu),sigma);
%                             case 'pfw', stat(neurons) = zscore(stat(neurons));
%                         end
%                 end
%             else
                stat(neurons) = zscore(stat(neurons));
%            end
                           
            switch stabilityType
                case 'time'           
                    %Get correlation coefficients and p-values. 
                    corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons);
                    %tuningStatus = TCRemap(mds(ssns(s)),mds(ssns(s+1)));
                    
                case 'place'
                    %Get the correlation coefficients and p-values.
                    corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons); 
                    %tuningStatus = PCRemap(mds(ssns(s)),mds(ssns(s+1))); 
            end
            
            %[~,stblcrit] = fdr_bh(corrStats(~isnan(corrStats(:,2)),2));
            stblcrit = .01/length(neurons);

            %Stable cells based on correlation.
            stable{a}{s} = intersect(find(corrs(:,2) < stblcrit),neurons);
            unstable{a}{s} = intersect(find(corrs(:,2) >= stblcrit | isnan(corrs(:,2))),neurons);
%             stable{a}{s} = intersect(find(tuningStatus(:,2)==1),neurons);
%             unstable{a}{s} = intersect(find(tuningStatus(:,2)<1),neurons);
            
            
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
        
        