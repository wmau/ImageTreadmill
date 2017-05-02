    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = MD(292:309); 
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
    
    [lags,pctMatches] = deal([]);
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        nSessions = length(ssns);
        
        mapMD = getMapMD(fulldataset(ssns)); 
        
        for s1=1:nSessions-1
            cd(fulldataset(ssns(s1)).Location);
            load('FinalOutput.mat','NumNeurons');
            neurons = 1:NumNeurons;
%             neurons = getTimeCells(fulldataset(ssns(s1)));
%             NumNeurons = length(neurons);
            
            for s2=s1+1:nSessions
                matches = msMatchCells(mapMD,[fulldataset(ssns(s1)) fulldataset(ssns(s2))],...
                    neurons,true);
                
                pctMatches = [pctMatches size(matches,1)/NumNeurons];
                lags = [lags s2-s1];
            end
        end
    end
    
    m = accumarray(lags',pctMatches',[],@mean);
    sem = accumarray(lags',pctMatches',[],@std)./accumarray(lags',pctMatches',[],@length);
    figure('Position',[520   518   334   280]);
    errorbar(1:4,m,sem)
    xlim([0,5])
    xlabel('Lag');
    ylabel('Proportion Overlap');
    ylabel('Proportion Active on Both Days');
        