function [flatDur,flatStats,m,sem,tukey] = InformationOverDaysStable(mds,stabilityType,statType,c)
%
%
%

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals);
    statType = lower(statType);
    stabilityType = lower(stabilityType); 
    
%% Determine what type of neuron to analyze.
    switch stabilityType
        case 'time'
            switch statType
                case 'ti',cellGet = 'timecells'; normCell = 'timecells';
                case 'si',cellGet = 'timecells'; normCell = 'placecells';
            end
        case 'place'
            switch statType
                case 'ti',cellGet = 'placecells'; normCell = 'timecells';
                case 'si',cellGet = 'placecells'; normCell = 'placecells';
            end
    end

%% Get a mapping matrix for neurons without recounting.
    [MAP,stabilityDuration,stats,stillthere] = deal(cell(nAnimals,1));
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{mds.Animal}));
        
        %Load registration matrix.
        mapMD = getMapMD(mds(ssns));
        cd(mapMD.Location);
        load('batch_session_map.mat');
        MAP{a} = batch_session_map.map(:,2:end);
        
        mapCols = zeros(length(ssns)-1,1);
        allrows = [];
        for s=1:length(ssns)-1
            cd(mds(ssns(s)).Location);
            
            neurons = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
            
            [~,mapRows,mapCols(s)] = msMatchCells(mapMD,mds(ssns(s)),neurons,false);
            allrows = [allrows; mapRows];
        end
        
        allrows = unique(allrows);
        
        %Matrix containing all our cells for this animal.
        MAP{a} = MAP{a}(allrows,mapCols);
        [stabilityDuration{a},stats{a},stillthere{a}] = deal(nan(size(MAP{a})));
 
%% Determine # days stable. 
        for s=1:length(ssns)-1
            cd(mds(ssns(s)).Location);
 
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI','Ispk','Isec');
                    stat = MI;                
                    
                case 'si'
                    load('SpatialInfo.mat','MI','Ispk','Isec');
                    stat = MI;
                
%                 case 'pk'
%                     [~,stat] = getTimePeak(mds(ssns(s)));
% 
%                     if strcmp(stabilityType,'time')
%                         neurons = TCs;
%                     end           
                
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [~,f] = size(PSAbool);
                    stat = sum(PSAbool,2)./f; 
                
            end
            
            %Neurons we are looking at.
            neurons = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
            
            %Neurons we are comparing to.
            norm = AcquireTimePlaceCells(mds(ssns(s)),normCell);
            
            %Get cells to normalize to. 
            norm = EliminateUncertainMatches(mds(ssns),norm);
            neurons = EliminateUncertainMatches(mds(ssns),neurons);
            
            %Normalize to a specified population.
%             m = mean(stat(norm)); 
%             sd = std(stat(norm));
%             stat(neurons) = bsxfun(@rdivide, bsxfun(@minus, stat(neurons), m), sd);
             
            %Normalize to themselves.
            %stat(neurons) = zscore(stat(neurons));
            
            %Normalize to entire population.
            stat = zscore(stat);
 
            %Neurons whose stability we are examining.
            %toCorr = intersect(neurons,MAP{a}(:,s));            
           
            [good,idx] = ismember(neurons,MAP{a}(:,s));
            idx(idx==0) = [];
            stats{a}(idx,s) = stat(neurons(good)); 
            
            %Set days stable for valid neurons to 0.
            if length(neurons) > 1
                stabilityDuration{a}(idx,s) = 0;
            else 
                stabilityDuration{a}(idx,s) = nan;
            end
            stillthere{a}(idx,s) = true;
%% 
            for ss=s+1:length(ssns)
                %corrMe = EliminateUncertainMatches([mds(ssns(s)),mds(ssns(ss))],toCorr);
                
                %Correlations.
                switch stabilityType
                    case 'time'
                        corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(ss)),neurons);               
                    case 'place'
                        corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(ss)),neurons);
                end
                
                %Stability criterion.
                stblcrit = .01/length(neurons);

                %Find stable neurons. 
                stableNeuron = find(corrs(:,2)<stblcrit);
                unstableNeuron = find(corrs(:,2)>stblcrit | isnan(corrs(:,2)));
                [~,idx] = ismember(stableNeuron,MAP{a}(:,s));
                idx(idx==0) = [];
                [~,uidx] = ismember(unstableNeuron,MAP{a}(:,s));
                uidx(uidx==0) = [];
                
                %Add a day to stable neurons. 
                ind = idx(idx & stillthere{a}(idx,s));
                stabilityDuration{a}(ind,s) = stabilityDuration{a}(ind,s)+1;
                stillthere{a}(uidx,s) = false;
                
            end
            

        end
     
        
    end
    
    flatDur = cell2mat(cellfun(@(x) x(:),stabilityDuration,'unif',false));
    flatStats = cell2mat(cellfun(@(x) x(:),stats,'unif',false));
    flatDur(isnan(flatDur)) = [];
    flatStats(isnan(flatStats)) = [];
    
    [flatDur,order] = sort(flatDur);
    flatStats = flatStats(order); 
    
%     edges = 0:4;
%     N = histc(flatDur,edges); 
%     badbin = edges(N<10);
%     bad = ismember(flatDur,badbin);
%     
%     flatDur(bad) = [];
%     flatStats(bad) = [];
     
    m = accumarray(flatDur+1,flatStats,[],@mean);
    sem = accumarray(flatDur+1,flatStats,[],@(x) std(x)/sqrt(length(x)));
    errorbar(unique(flatDur),m,sem,'linewidth',2,'color',c,'capsize',0,'marker','o');
    xlabel('Days Stable');
    ylabel('Mutual Information [z-scored bits]');
    set(gca,'tickdir','out');
    
end