%% Set up.
    clear; 
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309);      

    %Colors.
    teal = [0 .5 .5];
    purple = [.58 .44 .86];
    
%% 
    quartile = 'upper';
    cellType = 'place';
    infoType = 'si';
    
    infoDeltas(fulldataset,cellType,infoType);
    
function infoDeltas(mds,cellType,infoType)
%
%
%

%%
    animals = unique({mds.Animal}); 
    nAnimals = length(animals); 
    infoType = lower(infoType); 
    
    %Preallocate. 
    [sessionPairInfos,infoDelta] = deal(cell(nAnimals,1)); 
    [allInfoDeltas,originalInfos] = deal([]);
    %For each animal, get all the session pairs. 
    for a=1:nAnimals 
        ssns = find(strcmp(animals{a},{mds.Animal}));
        nSessions = length(ssns)-1;
        
        %Preallocate.
        [sessionPairInfos{a},infoDelta{a}] = deal(cell(nSessions,1));
        
        for s=1:nSessions
            %Reference and comparator session. 
            sessionPair = [mds(ssns(s)) mds(ssns(s+1))];
            
            %Get all non-coding neurons and eliminate those that didn't map
            %to the second session. 
            cd(mds(ssns(s)).Location); 
            load('FinalOutput.mat','NumNeurons'); 
            neurons = 1:NumNeurons; 
            neurons = EliminateUncertainMatches(sessionPair,neurons); 
            
            %Grab the mutual information scores for each neuron across a
            %day. 
            sessionPairInfos{a}{s} = msStats(sessionPair,infoType,neurons);
                      
            %Get the coding cells (time or place) and get the set
            %difference between them and all the cells. This gets you the
            %non-coding cells. 
            switch cellType
                case 'time', codingCells = getTimeCells(mds(ssns(s))); 
                case 'place', codingCells = getPlaceCells(mds(ssns(s)),.01);
            end
            nonCodingCells = ismember(1:NumNeurons,setdiff(neurons,codingCells))'; 
            
            %Get changes in information. 
            infoDelta{a}{s} = diff(sessionPairInfos{a}{s}(nonCodingCells,:),[],2);
            
            originalInfos = [originalInfos; sessionPairInfos{a}{s}(nonCodingCells,1)];
        end
        
         
        allInfoDeltas = [allInfoDeltas; cell2mat(infoDelta{a})];
    end
    scatter(originalInfos,,allInfoDeltas); 
    xlabel('Day 0 Info'); 
    ylabel('Delta Info'); 
end