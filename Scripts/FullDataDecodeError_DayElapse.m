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

daylag0_allErrors = nan(40,nSessions); 
daylag0_allErrorsShuffle = nan(40,B*nSessions);
sIndex = 1; c = 1;
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)
        %Get the decode error for each run. 
        [Mdl,~,testX] = TimeDecoder(fulldataset(ssns(s))); 
        decodedTime = PredictTime(Mdl,testX,'plotit',false);
        decodeError = TimeDecodeError(decodedTime); 
        
        %Take the mean across runs. 
        daylag0_allErrors(:,c) = mean(decodeError,2);
        
        %For each session, run B shuffle tests. 
        for i=1:B
            [badMdl,~,badTestX] = TimeDecoder(fulldataset(ssns(s)),'shuffle',true);
            shuffleDecode = PredictTime(badMdl,badTestX,'plotit',false);
            decodeErrorShuffle = TimeDecodeError(shuffleDecode); 

            daylag0_allErrorsShuffle(:,sIndex) = mean(decodeErrorShuffle,2); 
            sIndex = sIndex+1; 
        end
        c = c+1; 
    end
end

%Mean error from trial decodes, then the mean of that = the mean of the
%session. 
mean_error = cell2mat(cellfun(@(x) mean(x(:)),allErrors,'unif',0)); 
sem_error = cell2mat(cellfun(@(x) std(x(:))/sqrt(length(x(:))),allErrors,'unif',0)); 

% Results from day lag 0.
daylag0_errors = mean(daylag0_allErrors);
daylag0_error_mean = mean(daylag0_errors);
daylag0_error_sem = standarderror(daylag0_errors); 

mean_error = [daylag0_error_mean; mean_error];
sem_error = [daylag0_error_sem; sem_error];

%Same with the shuffles. 
shuffleM = cell2mat(cellfun(@(x) mean(x(:)),allErrorsShuffle,'unif',0));
shuffleSEM = cell2mat(cellfun(@(x) std(x(:))/sqrt(length(x(:))),...
    allErrorsShuffle,'unif',0));

daylag0_shuffleErrors = mean(decodeErrorShuffle(:));
daylag0_errorSEM = standarderror(decodeErrorShuffle(:));

shuffleM = [daylag0_shuffleErrors; shuffleM];
shuffleSEM = [daylag0_errorSEM; shuffleSEM];

%Do ANOVA of errors across days. 
[X,dayLags,condition] = deal([]); 
nComparisons = size(allErrors,1); 
for i=1:nComparisons
    X = [X, mean(allErrors{i})];
    
    dayLags = [dayLags, i.*ones(1,size(allErrors{i},2))];
    
    condition = [condition, ones(1,size(allErrors{i},2))];
end
X = [daylag0_errors, X];
dayLags = [zeros(1,length(daylag0_errors)), dayLags];
condition = [ones(1,length(daylag0_errors)), condition];

for i=1:nComparisons
    X = [X, mean(allErrorsShuffle{i})];
    
    dayLags = [dayLags, i.*ones(1,size(allErrorsShuffle{i},2))];
    
    condition = [condition, zeros(1,size(allErrorsShuffle{i},2))];
end
X = [mean(daylag0_allErrorsShuffle), X];
dayLags = [zeros(1,size(daylag0_allErrorsShuffle,2)), dayLags];
condition = [zeros(1,size(daylag0_allErrorsShuffle,2)), condition];

% [X_real,dayLags_real,condition_real] = ...
%     deal([]);
% nComparisons = size(allErrors,1); 
% for i=1:nComparisons
%     X_real = [X_real, mean(allErrors{i})];
%     
%     dayLags_real = [dayLags_real, i.*ones(1,size(allErrors{i},2))];
% 
%     condition_real = [condition_real, ones(1,size(allErrors{i},2))];
% end

%Also hard-coded, data from FullDataDecodeError, lines 45 and 46.
% X_real = [2.4379 2.1100 2.3249 2.1744 1.1233 1.2356 1.1140 1.3191 1.1282,...
%     1.5919 1.3403 2.1808 1.8733 3.1078 5.0681 2.4058, X_real];
% dayLags_real = [zeros(1,16), dayLags_real];
[p,tab,stats] = anovan(X,{dayLags, condition});

figure('Position',[680   428   320   550]); hold on;
errorbar([0:3],mean_error,sem_error,'k','linewidth',2);
errorbar([0:3],shuffleM,shuffleSEM,'r','linewidth',2); 
xlim([-.5 3.5]);
set(gca,'tickdir','out','xtick',(0:3));
xlabel('Day lag'); 
ylabel('Decode error (s)');

