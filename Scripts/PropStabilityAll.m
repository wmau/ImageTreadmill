%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = [MD(292:303) MD(305:308)]; 
   
    cellType = 'time';
    switch cellType
        case 'time', c = [0 .5 .5];
        case 'place', c = [.58 .44 .86];
    end
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
       
    [map,onlyIncoming,onlyOutgoing,bothIandO] = deal(cell(nAnimals,1));
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        [pct(a),stability(a),map{a}] = PropStability2(fulldataset(ssns),cellType);
        
        onlyIncoming{a} = setdiff(stability(a).Incoming,stability(a).Outgoing);
        onlyOutgoing{a} = setdiff(stability(a).Outgoing,stability(a).Incoming);
        bothIandO{a} = intersect(stability(a).Outgoing,stability(a).Incoming);
    end
    
    sPct = round(mean([pct.Stable])*100,1);
    iPct = round(mean(cellfun('length',onlyIncoming) ./ cellfun('length',map)).*100,1);
    oPct = round(mean(cellfun('length',onlyOutgoing) ./ cellfun('length',map)).*100,1);
    ioPct = round(mean(cellfun('length',bothIandO) ./ cellfun('length',map)).*100,1);
    
    
    figure;
    p = pie([sPct iPct ioPct oPct],[1 0 0 0],...
        {   ['Stable (',num2str(sPct),'%)'],...
            ['Entering (',num2str(iPct),'%)'],...
            ['Transient (',num2str(ioPct),'%)'],...
            ['Exiting (',num2str(oPct),'%)']});
    p(1).FaceColor = c;
    p(3).FaceColor = c;
    p(3).FaceAlpha = .3;
    p(5).FaceColor = [.5 .5 .5];
    p(7).FaceColor = [0 0 0];
