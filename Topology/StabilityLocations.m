function [D,stabilityStatus] = StabilityLocations(mds,stabilityType)
%
%
%

%%
    animals = unique({mds.Animal});
    nAnimals = length(animals); 
    stabilityType = lower(stabilityType);
    sf = 1.1;
    
%% What type of neuron to analyze. 
    switch stabilityType
        case 'time',cellGet = 'timecells'; 
        case 'place',cellGet = 'placecells';
    end
    
%% 
    [D,stabilityStatus.stable,stabilityStatus.unstable] = deal(cell(nAnimals,1));
    for a=1:nAnimals
        %Get all the sessions for this animal.
        ssns = find(strcmp(animals{a},{mds.Animal})); 
        nSessions = length(ssns)-1;
       
        [D{a},stabilityStatus.stable{a},stabilityStatus.unstable{a}] = deal(cell(nSessions,1));
        for s=1:nSessions
            cd(mds(ssns(s)).Location);
            load('FinalOutput.mat','NumNeurons');
            D{a}{s} = nan(NumNeurons);
            
            %Get neurons and their centroids. 
            neurons = AcquireTimePlaceCells(mds(ssns(s)),cellGet);
            centroids = getNeuronCentroids(mds(ssns(s)));
            
            %Calculate distance.
            for n1=1:NumNeurons
                %Centroid for cell 1.
                x1 = centroids(n1,1);
                y1 = centroids(n1,2);
                
                for n2=n1+1:NumNeurons
                    %Centroid for cell 2.
                    x2 = centroids(n2,1);
                    y2 = centroids(n2,2); 
                    
                    %Anatomical distance. 
                    D{a}{s}(n1,n2) = sqrt((x2-x1)^2 + (y2-y1)^2)*sf; 
                end
            end
                
            %Determine stability. 
            switch stabilityType
                case 'time'
                    corrs = CorrTrdmllTrace(mds(ssns(s)),mds(ssns(s+1)),neurons);
                case 'place'
                    corrs = CorrPlaceFields(mds(ssns(s)),mds(ssns(s+1)),neurons);
            end
            
            stblcrit = .01/length(neurons);
            
            good = find(corrs(:,2) < stblcrit);
            bad = find(corrs(:,2) > stblcrit | isnan(corrs(:,2)));
            
            stabilityStatus.stable{a}{s} = intersect(good,neurons);
            stabilityStatus.unstable{a}{s} = intersect(bad,neurons); 
         
        end
    end
    
end