%%
%Posterior probability density for all sessions and all trials averaged
%together. 

%%
clear;
loadMD;

fulldataset = [MD(292:303) MD(305:308)];
nSessions = length(fulldataset); 

animals = unique({fulldataset.Animal});
nAnimals = length(animals); 

postProbs = [];
for a=1:nAnimals
    ssns = find(strcmp(animals{a},{fulldataset.Animal}));
    
    for s=1:length(ssns)
        %Get the decode error for each run. 
        [Mdl,~,testX,testLaps] = TimeDecoder(fulldataset(ssns(s))); 
        [decodedTime,temp] = PredictTime(Mdl,testX,'plotit',false);
        
        postProbs = cat(3,postProbs,temp);
    end
    
end

TimeDecoderPosteriorProbabilityPlot(postProbs); 