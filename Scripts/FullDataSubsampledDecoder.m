clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

B = 10;
pctCells = 0.05:0.1:1;
n = length(pctCells); 
performance = zeros(nSessions,n); 
c = 1;
[control] = deal([]);
for i=1:n
    shufflePerformances = [];
    for a=1:nAnimals
        ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
        
        for s=1:length(ssns)
            TCs = getTimeCells(fulldataset(ssns(s)));
            nTCs = ceil(length(TCs)*pctCells(i));
            sampleTCs = randsample(TCs,nTCs); 
            
            [Mdl,~,testX] = TimeDecoder(fulldataset(ssns(s)),'neurons',...
                sampleTCs); 
            [decodedTime,postProbs] = PredictTime(Mdl,testX,'plotit',false); 
            decodeError = TimeDecodeError(decodedTime);
            
            performance(c,i) = mean(decodeError(:));
            
            
            for j=1:B
                [badMdl,~,badTestX] = TimeDecoder(fulldataset(ssns(s)),...
                    'shuffle',true,'neurons',sampleTCs);
                shuffleDecode = PredictTime(badMdl,badTestX,'plotit',false);
                decodeErrorShuffle = TimeDecodeError(shuffleDecode);
               
                shufflePerformances = [shufflePerformances; mean(decodeErrorShuffle(:))];
            end
            c = c+1;
        end    
    end
    
    control = [control shufflePerformances];
    c = 1;
end


f = figure;
errorbar(pctCells,mean(performance),standarderror(performance));
xlabel('Proportion of cells in classifier');
ylabel('Decode error (s)');
make_plot_pretty(gca);