function [STATS,stabilityStatus] = StabilityStats(mds,stabilityType,statType)
%
%
%

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals); 
    statType = lower(statType);
    
%%
    %What type of neuron to analyze. 
    switch statType
        case 'ti'
            switch stabilityType
                case 'time',cellGet = 'timecells';
                case 'place',cellGet = 'dual';
            end
        case 'si'
            switch stabilityType
                case 'time',cellGet = 'dual';
                case 'place',cellGet = 'placecells';
            end
        case {'fr','fluor'}
            switch stabilityType
                case 'time',cellGet = 'timecells';
                case 'place',cellGet = 'placecells';
            end
        case 'pk'
            cellGet = 'timecells'; 
    end

%% Only look at each cell once. 
    [map,neurons] = deal(cell(nAnimals,1));
    for a=1:nAnimals 
        %Get all the sessions for this animal except the last one. 
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        ssns(end) = []; 
        nSessions = length(ssns);
        
        neurons{a} = cell(nSessions,1);
        for s=1:nSessions
            neurons{a}{s} = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
        end
        
        %Make a matrix containing each cell of interest, without
        %recounting.
        map{a} = msMatchMultiSessionCells(mds(ssns),neurons{a});
    end

%%
    [STATS.stable,STATS.unstable] = deal(cell(1,nAnimals));
    [nNeurons.stable] = deal(zeros(1,nAnimals)); 
    
    for a=1:nAnimals       
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions = length(ssns)-1;
       
        for s=1:nSessions
            cd(mds(ssns(s)).Location);
            
            N = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
            n = intersect(N,map{a}(:,s));
            
            switch statType
                %Temporal information.
                case 'ti'
                    load('TemporalInfo.mat','MI','Ispk','Isec');
                    stat = MI; 
                    
                %Spatial information.
                case 'si'
                    load('SpatialInfo.mat','MI','Ispk','Isec');
                    stat = MI;
                    
                %Firing rate (really Ca event rate).
                case 'fr'
                    load('Pos_align.mat','PSAbool');
                    [~,f] = size(PSAbool);
    %                 d = diff([zeros(n,1) PSAbool],1,2);
    %                 d(d<0) = 0;
                    stat = sum(PSAbool,2)./f; 
                    stat = (stat-min(stat))./range(stat);
                    
                %Mean fluorescence.
                case 'fluor'
                    load('Pos_align.mat','DFDTtrace');
                    stat = mean(DFDTtrace,2);
                                    
                %Time peak. 
                case 'pk'
                    [~,stat] = getTimePeak(mds(ssns(s)));    
            end
         
            switch stabilityType
                case 'time'           
                    %Get correlation coefficients and p-values. 
                    corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),n);

                case 'place'
                    %Get the correlation coefficients and p-values.
                    corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),n);
            end
            
            %Normalize.
            if ~any(strcmp(statType,{'fr','fluor'}))
                stat(N) = zscore(stat(N));
            end
            
            stblcrit = .05/length(N);
           
            good = find(corrs(:,2) < stblcrit);
            bad = find(corrs(:,2) > stblcrit | isnan(corrs(:,2))); 
            
            validStable = intersect(good,n);
            validUnstable = intersect(bad,n); 
            
            stabilityStatus.stable{a}{s} = validStable;
            stabilityStatus.unstable{a}{s} = validUnstable;
            
            STATS.stable{a}{s} = stat(validStable); 
            STATS.unstable{a}{s} = stat(validUnstable);            
        end
                    
        
    end
end