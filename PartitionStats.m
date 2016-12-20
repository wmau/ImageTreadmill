function [STATS,nNeurons,stable,unstable] = PartitionStats(mds,stabilityCriterion,statType)
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

%% Compile
    STATS.stable = cell(1,nAnimals);
    STATS.unstable = cell(1,nAnimals);
    nNeurons.stable = zeros(1,nAnimals); 
    nNeurons.unstable = zeros(1,nAnimals);
    stable = cell(1,nAnimals);
    unstable = cell(1,nAnimals);
    for a = 1:nAnimals       
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        
        nStable = 0; 
        nUnstable = 0;
        STATS.stable{a} = [];
        STATS.unstable{a} = [];
        stable{a} = cell(1,length(ssns)-1);
        unstable{a} = cell(1,length(ssns)-1);
        for s = 1:length(ssns)-1
            cd(mds(ssns(s)).Location);
            
            load('TimeCells.mat','TimeCells');
            load('TemporalInfo.mat','sig'); 
            load('Placefields.mat','pval');
            load('PlacefieldStats.mat','PFnHits','PFpcthits','bestPF');
            
            PCcrit = .01;
            stblcrit = .05;
            
            %Get all time cells with a viable place field. 
            idx = sub2ind(size(PFnHits), 1:size(PFnHits,1), bestPF');
            
            if strcmp(statType,'TI')
                load('TemporalInfo.mat','MI','Ispk','Isec');
                stat = MI; 
                noi = intersect(find(sig),TimeCells);
            elseif strcmp(statType,'SI')
                load('SpatialInfo.mat','MI','Ispk','Isec');
                stat = MI;
                noi = find(pval<PCcrit & PFnHits(idx) > 4);
            elseif strcmp(statType,'FR')
%                 load('Pos_align.mat','FT');
%                 [n,f] = size(FT);
%                 d = diff([zeros(n,1) FT],1,2);
%                 d(d<0) = 0;
%                 stat = sum(d,2)./f; 
                load(fullfile(pwd,'TimeCells.mat'),'TodayTreadmillLog','T');
                load(fullfile(pwd,'Pos_align.mat'),'PSAbool');
                inds = TrimTrdmllInds(TodayTreadmillLog,T);
                n = size(FT,1);
                rasters = cell(1,n);
                stat = zeros(n,1);
                for nn=1:n
                    rasters{nn} = buildRaster(inds,FT,nn,'onsets',false);
                    stat(nn) = sum(rasters{nn}(:))./numel(rasters{nn});
                end
            end
                    
            if strcmp(stabilityCriterion,'time')
                noi = intersect(noi,intersect(find(sig),TimeCells));
                
                %Get correlation coefficients and p-values. 
                corrStats = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),noi);
                tuningStatus = TCRemap(mds(ssns(s)),mds(ssns(s+1)));

                %Stable time cells based on correlation and non-shifting time
                %field.
                stable{a}{s} = intersect(find(corrStats(:,2) < stblcrit & tuningStatus(:,2)==1),noi);
                unstable{a}{s} = intersect(find(corrStats(:,2) >= stblcrit | tuningStatus(:,2)~=1),noi);
         
            elseif strcmp(stabilityCriterion,'place')
                noi = intersect(noi,find(pval<PCcrit & PFnHits(idx) > 4));
                
                %Get the correlation coefficients and p-values.
                corrStats = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),noi);

                %Stable place cells based on correlation p-value.
                stable{a}{s} = intersect(find(corrStats(:,2) < stblcrit),noi); 
                unstable{a}{s} = intersect(find(corrStats(:,2) >= stblcrit | isnan(corrStats(:,2))),noi);     
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
        
        