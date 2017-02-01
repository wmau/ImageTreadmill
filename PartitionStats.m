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
    for a = 1:nAnimals       
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        
        [nStable,nUnstable] = deal(0); 
        [STATS.stable{a},STATS.unstable{a}] = deal([]);
        [stable{a},unstable{a}] = deal(cell(1,length(ssns)-1));
        for s = 1:length(ssns)-1
            cd(mds(ssns(s)).Location);
            
            PCcrit = .01;
            
            %Get place and time cells.
            PCs = getPlaceCells(mds(ssns(s)),PCcrit);
            TCs = getTimeCells(mds(ssns(s)));
                        
            if strcmp(statType,'ti')
                load('TemporalInfo.mat','MI','Ispk','Isec');
                stat = MI;                
                
                if strcmp(stabilityType,'time')
                    neurons = TCs;
                elseif strcmp(stabilityType,'place')
                    neurons = intersect(TCs,PCs);
                end
            elseif strcmp(statType,'si')
                load('SpatialInfo.mat','MI','Ispk','Isec');
                stat = MI;
                
                if strcmp(stabilityType,'time')
                    neurons = intersect(TCs,PCs);
                elseif strcmp(stabilityType,'place')
                    neurons = PCs;
                end
                
            elseif strcmp(statType,'pk')
                [~,stat] = getTimePeak(mds(ssns(s)));
                
                if strcmp(stabilityType,'time')
                    neurons = TCs;
                end           
                
            elseif strcmp(statType,'fr')
                load('Pos_align.mat','PSAbool');
                [n,f] = size(PSAbool);
%                 d = diff([zeros(n,1) PSAbool],1,2);
%                 d(d<0) = 0;
                stat = sum(PSAbool,2)./f; 
                
                if strcmp(stabilityType,'time')
                    neurons = TCs;
                elseif strcmp(stabilityType,'place')
                    neurons = PCs;
                end

%                 load('Pos_align.mat','DFDTtrace');
%                 stat = mean(DFDTtrace,2);

%                 load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','T');
%                 load(fullfile(pwd,'Pos_align.mat'),'PSAbool');
%                 inds = TrimTrdmllInds(TodayTreadmillLog,T);
%                 n = size(PSAbool,1);
%                 rasters = cell(1,n);
%                 stat = zeros(n,1);
%                 for nn=1:n
%                     rasters{nn} = buildRaster(inds,PSAbool,nn,'onsets',false);
%                     stat(nn) = sum(rasters{nn}(:))./numel(rasters{nn});
%                 end                       
            end
            
            %stat = (stat-min(stat))./range(stat);
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
    end
    
end
        
        