clear; 
loadMD; 

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 100;
[pCorrectIterations,pCorrectStable,pCorrectUnstable] = ...
    deal(cell(1,nAnimals));
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal})); 
    nSessions = length(ssns);
    
    [~,stability,map,everStable] = PropStability2(fulldataset(ssns),'time');
    activeOnAllDays = all(map > 0 & ~isnan(map),2);
    stableCellIndicesOnMap = find(everStable & activeOnAllDays);
    %stableCellIndicesOnMap = intersect(find(activeOnAllDays),stability.Stable);
    stableCells = cell(nSessions,1);
    for s=1:nSessions
        stableCells{s} = map(stableCellIndicesOnMap,s);
    end
    
    nStable = length(stableCellIndicesOnMap);
    nUnstable = length(setdiff(find(activeOnAllDays),find(everStable)));
    
    disp(['# of stable cells:',num2str(nStable)]);
    disp(['# of unstable cells:',num2str(nUnstable)]);
    
    [pCorrectIterations{a},pCorrectStable{a},pCorrectUnstable{a}] = ...
        deal(nan(1,B)); 
    for i=1:B
%         %Get decode error for each run. 
%         [Mdl,~,testX,~,testY] = DayDecoder(fulldataset(ssns)); 
%         [decodedDay,~] = PredictDay(Mdl,testX); 
%         pCorrectIterations{a}(i) = DayDecodeError(decodedDay,testY); 
        
        %Get accuracy of decoder trained without stable cells. 
        [Mdl,~,testX,~,testY] = DayDecoder(fulldataset(ssns),'neurons',...
            stableCells);
        [decodedDay,~] = PredictDay(Mdl,testX); 
        pCorrectStable{a}(i) = DayDecodeError(decodedDay,testY); 
        
        %Get accuracy of decoder trained without random cells. 
        randomCells = cell(nSessions,1);
        %randomSample = randsample(find(activeOnAllDays),...
            %length(unstableCellIndicesOnMap));
        randomSample = randsample(setdiff(find(activeOnAllDays),find(~everStable)),...
            nStable);
        for s=1:nSessions
            randomCells{s} = map(randomSample,s);
        end
        [Mdl,~,testX,~,testY] = DayDecoder(fulldataset(ssns),'neurons',...
            randomCells);
        [decodedDay,~] = PredictDay(Mdl,testX); 
        pCorrectUnstable{a}(i) = DayDecodeError(decodedDay,testY); 
    end
    
end

omitStable = cell2mat(pCorrectStable); 
omitUnstable = cell2mat(pCorrectUnstable); 

m = [mean(omitUnstable) mean(omitStable)];
sem = [standarderror(omitUnstable) standarderror(omitStable)];
h = figure('position',[680   380   245   598]);
hold on;
for a=1:nAnimals
    mThisAnimal = [mean(pCorrectUnstable{a}) mean(pCorrectStable{a})];
    semThisAnimal = [standarderror(pCorrectUnstable{a}) standarderror(pCorrectStable{a})];
    errorbar([1,2],mThisAnimal,semThisAnimal,'color',[.7 .7 .7]);
end
errorbar([1,2],m,sem,'color','k','linewidth',4);
xlim([.5 2.5]); 
make_plot_pretty(h); 
set(gca,'tickdir','out','xtick',[1,2],...
    'xticklabel',{'Stable cells','Unstable cells'});
ylabel('% trials correctly decoded');
