clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 50;
[allErrors,allErrorsShuffle] = deal(cell(3,1)); 
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    %Get all combinations. 
    combinations = combnk(ssns,2); 
    nCombs = size(combinations,1); 
    
    %For each session combination (within an animal), decode across days. 
    for c = 1:nCombs
        %Get sessions.
        training = combinations(c,1); 
        test = combinations(c,2); 
        dayLag = test-training; 
        
        %Get stable and unstable cells. 
        trainingTCs = getTimeCells(fulldataset(training));
        corrStats = CorrTrdmllTrace(fulldataset(training),fulldataset(test),...
            trainingTCs);
        pCrit = 0.01;
        %[~,pCrit] = fdr_bh(corrStats(~isnan(corrStats(:,1)),2));
        stableCells = find(corrStats(:,2) < pCrit);
        unstableCells = find(corrStats(:,2) >= pCrit);
        
        if length(unstableCells) > length(stableCells)
            unstableCells = randsample(unstableCells,length(stableCells));
        end
        
        %Run decoder. 
        [decodedTime,~,Mdl] = ElapsedDayTimeDecoder(fulldataset(training),...
            fulldataset(test),'plotit',false,'neurons',unstableCells);
        decodeError = TimeDecodeError(decodedTime); 
        allErrors{dayLag} = [allErrors{dayLag} mean(decodeError,2)];
   
        for i=1:B
            sample = randsample(stableCells,length(unstableCells));
            randomUnstableDecodedTime = ElapsedDayTimeDecoder(fulldataset(training),...
                fulldataset(test),'plotit',false,'Mdl',Mdl,'neurons',sample); 
            decodeErrorShuffle = TimeDecodeError(randomUnstableDecodedTime); 
            
            allErrorsShuffle{dayLag} = [allErrorsShuffle{dayLag} mean(decodeErrorShuffle,2)];
        end
    end
end

%Mean error from trial decodes, then the mean of that = the mean of the
%session. 
m = cell2mat(cellfun(@(x) mean(x(:)),allErrors,'unif',0)); 
sem = cell2mat(cellfun(@(x) std(x(:))/sqrt(length(x(:))),allErrors,'unif',0)); 

%Same with the shuffles. 
shuffleM = cell2mat(cellfun(@(x) mean(x(:)),allErrorsShuffle,'unif',0));
shuffleSEM = cell2mat(cellfun(@(x) std(x(:))/sqrt(length(x(:))),...
    allErrorsShuffle,'unif',0));

%Do ANOVA of error across days. 
X = [];
grps = [];
for i=1:3
    X = [X; allErrors{i}(:)];
    grps = [grps; i*ones(length(allErrors{i}(:)),1)];
end
[~,~,stats] = kruskalwallis(X,grps);
multcompare(stats);

figure('Position',[680   428   320   550]); hold on;
errorbar(m,sem,'k','linewidth',2);
errorbar(shuffleM,shuffleSEM,'r','linewidth',2); 
xlim([.5 3.5]);
set(gca,'tickdir','out','xtick',(1:3));
xlabel('Day lag'); 
ylabel('Decode error (s)');