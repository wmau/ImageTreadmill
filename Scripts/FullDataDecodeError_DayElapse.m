clear;
loadMD;

fulldataset = [MD(292:299) MD(300:303) MD(305:308)];
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
        training = combinations(c,1); 
        test = combinations(c,2); 
        dayLag = test-training; 
        [decodedTime,~,Mdl] = ElapsedDayTimeDecoder(fulldataset(training),...
            fulldataset(test),'plotit',false);
        decodeError = TimeDecodeError(decodedTime); 
        allErrors{dayLag} = [allErrors{dayLag} mean(decodeError,2)];
   
        for i=1:B
            shuffleDecode = ElapsedDayTimeDecoder(fulldataset(training),...
                fulldataset(test),'plotit',false,'shuffle',true); 
            decodeErrorShuffle = TimeDecodeError(shuffleDecode); 
            
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

%Do ANOVA of errors across days. 
[X_real,X_shuffle,dayLags_real,dayLags_shuffle,condition_real,condition_shuffle] = ...
    deal([]);
nComparisons = size(allErrors,1); 
for i=1:nComparisons
    X_real = [X_real, mean(allErrors{i})];
    X_shuffle = [X_shuffle, mean(allErrorsShuffle{i})];
    
    dayLags_real = [dayLags_real, i.*ones(1,size(allErrors{i},2))];
    dayLags_shuffle = [dayLags_shuffle, i.*ones(1,size(allErrorsShuffle{i},2))];

    condition_real = [condition_real, ones(1,size(allErrors{i},2))];
    condition_shuffle = [condition_shuffle, zeros(1,size(allErrorsShuffle{i},2))];
end
[p,tab,stats] = anovan([X_real X_shuffle],{[condition_real,condition_shuffle],...
    [dayLags_real dayLags_shuffle]});

figure('Position',[680   428   320   550]); hold on;
errorbar(m,sem,'k','linewidth',2);
errorbar(shuffleM,shuffleSEM,'r','linewidth',2); 
xlim([.5 3.5]);
set(gca,'tickdir','out','xtick',(1:3));
xlabel('Day lag'); 
ylabel('Decode error (s)');

