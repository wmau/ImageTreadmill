clear; 
loadMD; 

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 50;
percentCorrect = nan(1,nAnimals);
pCorrectShuffleIterations = nan(1,nAnimals*B);
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal})); 
    
    [pCorrectIterations] = deal(nan(1,B));
    for i=1:B
        %Get decode error for each runs. 
        [Mdl,~,testX,~,testY] = DayDecoder(fulldataset(ssns)); 
        [decodedDay,postProbs] = PredictDay(Mdl,testX); 
        pCorrectIterations(i) = DayDecodeError(decodedDay,testY); 
        
        %shuffleDecodedDay = PredictDay(Mdl,testX(:,randperm(size(testX,2))));
        shuffleDecodedDay = PredictDay(Mdl,testX);
        pCorrectShuffleIterations(i*a) = DayDecodeError(shuffleDecodedDay,testY(randperm(length(testY)))); 
    end
    
    percentCorrect(a) = mean(pCorrectIterations); 
end

scatterBox([percentCorrect,pCorrectShuffleIterations],...
    [zeros(size(percentCorrect)), ones(size(pCorrectShuffleIterations))],...
    'yLabel','% test trials correctly decoded','xLabels',...
    {'Real','Shuffle'},'boxColor',[0 0 0;1 0 0]);
set(gca,'fontsize',14);
