%% Set up.
    clear;
    
    %Load MD entries. 
    loadMD;
    %MD(292:295) = G45.
    %MD(296:299) = G48.
    %MD(300:304) = Bellatrix.
    %MD(305:309) = Polaris.
    fulldataset = [MD(292:299) MD(300:303) MD(305:308)]; 
    
    modality = 'time';
    
    animals = unique({fulldataset.Animal});
    nAnimals = length(animals);
       
    [map,onlyIncoming,onlyOutgoing,bothIandO] = deal(cell(nAnimals,1));
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
        
        [pct(a),stability(a),map{a}] = PropStability(fulldataset(ssns),modality);
        
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
            ['Incoming (',num2str(iPct),'%)'],...
            ['Both (',num2str(ioPct),'%)'],...
            ['Outgoing (',num2str(oPct),'%)']});
%     p(1).FaceColor = [0 .5 .5];
%     p(3).FaceColor = [0 .5 .5];
%     p(3).FaceAlpha = .5;