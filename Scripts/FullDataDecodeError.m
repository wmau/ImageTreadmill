clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 50;
[allErrors] = deal(nan(40,nSessions));
allErrorsShuffle = nan(40,B*nSessions); 
c = 1;
sIndex = 1; 
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)
        %Get the decode error for each run. 
        [Mdl,~,testX] = TimeDecoder(fulldataset(ssns(s))); 
        [decodedTime,postProbs] = PredictTime(Mdl,testX,'plotit',false);
        decodeError = TimeDecodeError(decodedTime); 
        
        %Take the mean across runs. 
        allErrors(:,c) = mean(decodeError,2);
        
        %For each session, run B shuffle tests. 
        for i=1:B
            [badMdl,~,badTestX] = TimeDecoder(fulldataset(ssns(s)),'shuffle',true);
            shuffleDecode = PredictTime(badMdl,badTestX,'plotit',false);
            decodeErrorShuffle = TimeDecodeError(shuffleDecode); 

            allErrorsShuffle(:,sIndex) = mean(decodeErrorShuffle,2); 
            sIndex = sIndex+1; 
        end
        c = c+1; 
    end
end

%% Plot error by time bin. 
%Time vector.
t = linspace(0,10,40);                     

%Mean decode error for time bin.
m = mean(allErrors,2)';                      
sem = (std(allErrors,[],2)/sqrt(nSessions))';  

%Mean decode error for shuffles. 
shuffleM = mean(allErrorsShuffle,2)';        
shuffleSEM = (std(allErrorsShuffle,[],2)/sqrt(nSessions*B))'; 

%Plot. 
empLine.col = {'k'};
shuffleLine.col = {'r'};

figure; hold on;
plot(t,allErrors,'color',[.6 .6 .6 .5]); 
mseb(t,m,sem,empLine,1);
mseb(t,shuffleM,shuffleSEM,shuffleLine,1); 
set(gca,'tickdir','out','xtick',[0 5 10],'ytick',[0 5 10]);
xlabel('Real time (s)');
ylabel('Decode error (s)');

m = median(allErrors,2);
m = repmat(m,1,B*nSessions);
p = sum(m > allErrorsShuffle,2) ./ (B*nSessions); 

sig = zeros(1,40);
sig(p < 0.05) = 1;
y = 10*ones(1,40);
y(~sig) = nan;
plot(t,y,'r');

%% Plot overall error. 
M = mean(allErrors); 

overallShuffleM = mean(allErrorsShuffle); 

scatterBox([M overallShuffleM],[zeros(size(M)) ones(size(overallShuffleM))],...
    'xLabels',{'Empirical','Shuffle'},'yLabel','Decode error (s)');

errorP = ranksum(

%% ANOVA on errors
[p,tab,stats] = anova1(allErrors);