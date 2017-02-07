function [flatDur,flatStats,m,sem,tukey] = InformationOverDaysStable(mds,stabilityType,statType)
%
%
%

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals);
    statType = lower(statType);
    stabilityType = lower(stabilityType); 
    PCcrit = .01;

%% 
    [MAP,stabilityDuration,stats] = deal(cell(nAnimals,1));
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
            
            TCs = getTimeCells(mds(ssns(s)));
            PCs = getPlaceCells(mds(ssns(s)),PCcrit);
            
            switch statType
                case 'ti'             
                    if strcmp(stabilityType,'time')
                        neurons = TCs;
                    elseif strcmp(stabilityType,'place')
                        neurons = intersect(TCs,PCs);
                    end  
                    
                case 'si'
                    if strcmp(stabilityType,'time')
                        neurons = intersect(TCs,PCs);
                    elseif strcmp(stabilityType,'place')
                        neurons = PCs;
                    end
                    
                case 'fr'
                    if strcmp(stabilityType,'time')
                        neurons = TCs;
                    elseif strcmp(stabilityType,'place')
                        neurons = PCs;
                    end
                    
            end
            
            [~,mapRows,mapCols(s)] = msMatchCells(mapMD,mds(ssns(s)),neurons,false);
            allrows = [allrows; mapRows];
        end
        
        allrows = unique(allrows);
        
        %Matrix containing all our cells for this animal.
        MAP{a} = MAP{a}(allrows,mapCols);
        [stabilityDuration{a},stats{a}] = deal(nan(size(MAP{a})));
 
%% 
        for s=1:length(ssns)-1
            cd(mds(ssns(s)).Location);
            
            TCs = getTimeCells(mds(ssns(s)));
            PCs = getPlaceCells(mds(ssns(s)),PCcrit);
            
            switch statType
                case 'ti'
                    load('TemporalInfo.mat','MI','Ispk','Isec');
                    stat = MI;                

                    if strcmp(stabilityType,'time')
                        neurons = TCs;
                    elseif strcmp(stabilityType,'place')
                        neurons = intersect(TCs,PCs);
                    end
                    
                case 'si'
                    load('SpatialInfo.mat','MI','Ispk','Isec');
                    stat = MI;

                    if strcmp(stabilityType,'time')
                        neurons = intersect(TCs,PCs);
                    elseif strcmp(stabilityType,'place')
                        neurons = PCs;
                    end
                
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
                
                    if strcmp(stabilityType,'time')
                        neurons = TCs;
                    elseif strcmp(stabilityType,'place')
                        neurons = PCs;
                    end
            end
            
            corrMe = intersect(neurons,MAP{a}(:,s));
            corrMe = EliminateUncertainMatches([mds(ssns(s)),mds(ssns(s+1))],corrMe);
            
            if strcmp(statType,'fr')
                stat(neurons) = (stat(neurons)-min(stat(neurons)))./range(stat(neurons));
            else
                stat(neurons) = zscore(stat(neurons));
            end
            
            [good,idx] = ismember(neurons,MAP{a}(:,s));
            idx(idx==0) = [];
            stats{a}(idx,s) = stat(neurons(good));  
            
            %Set days stable for valid neurons to 0.
            stabilityDuration{a}(idx,s) = 0;
            
%% 
            for ss=s+1:length(ssns)
                %Correlations.
                switch stabilityType
                    case 'time'
                        corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(ss)),corrMe);               
                    case 'place'
                        corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(ss)),corrMe);
                end
                
                %Stability criterion.
                stblcrit = .01/length(neurons);

                %Find stable neurons. 
                stableNeuron = find(corrs(:,2)<stblcrit);
                [~,idx] = ismember(stableNeuron,MAP{a}(:,s));
                idx(idx==0) = [];
                
                %Add a day to stable neurons. 
                stabilityDuration{a}(idx,s) = stabilityDuration{a}(idx,s)+1;
                
                
            end
            
            
        end
     
        
    end
    
    flatDur = cell2mat(cellfun(@(x) x(:),stabilityDuration,'unif',false));
    flatStats = cell2mat(cellfun(@(x) x(:),stats,'unif',false));
    flatDur(isnan(flatDur)) = [];
    flatStats(isnan(flatStats)) = [];
    %m = accumarray(flatDur+1,flatStats,[],@mean);
    
    [flatDur,order] = sort(flatDur);
    flatStats = flatStats(order);
    
    [p,tbl,anovastats] = kruskalwallis(flatStats,flatDur,'on');
    tukey = multcompare(anovastats,'display','on');
     
    m = accumarray(flatDur+1,flatStats,[],@mean);
    sem = accumarray(flatDur+1,flatStats,[],@(x) std(x)/sqrt(length(x)));
    errorbar(unique(flatDur),m,sem,'linewidth',2);
    xlabel('Days Stable');
    ylabel('Mutual Information [z-scored bits]');
    set(gca,'tickdir','out');
    
end