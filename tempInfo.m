function [I,sig] = tempInfo(MD)
%[I,sig] = tempInfo(MD)
%
%   Calculates the temporal information of all the neurons in a session. 
%
%   INPUT
%       MD: Session entry. 
%
%   OUPUTS
%       I: Temporal information.
%
%       sig: Significantly high temporal information, compared to permuted
%       data. 
%

%% Set up.
    currdir = pwd; 
    cd(MD.Location); 
    load(fullfile(pwd,'TimeCells.mat'),'ratebylap','TodayTreadmillLog','T');
    
    %Number of shuffles. 
    B = 1000;
    good = TodayTreadmillLog.complete & TodayTreadmillLog.delaysetting == T;    %Only look at complete laps at duration T. 
    ratebylap = ratebylap(good,:,:); 
    
%% Compute temporal information and permutation test. 
    [nLaps,nBins,nNeurons] = size(ratebylap); 
    I = zeros(nNeurons,1);              %Preallocate empirical temporal information.
    surrogate = zeros(nNeurons,B);      %Preallocate surrogate data. 
    p = ProgressBar(nNeurons);
    for n=1:nNeurons
        %Empirical temporal information. 
        I(n) = tempInfoOneNeuron(ratebylap(:,:,n));
        
        %Preallocate permuted ratebylap. 
        ratebylapShift = zeros(nLaps,nBins);
        for i=1:B
            %Card shuffle. 
            jitters = randsample(nBins,nLaps,true);
            
            %Make permuted ratebylap. 
            for l=1:nLaps
                ratebylapShift(l,:) = circshift(ratebylap(l,:,n),[0,jitters(l)]);
            end
            
            %Compute temporal information of permuted ratebylap. 
            surrogate(n,i) = tempInfoOneNeuron(ratebylapShift);
        end
    
        p.progress;
    end
    p.stop;
    
    %P-value of temporal information being significantly better than
    %shuffled data. 
    pval = sum(surrogate>repmat(I,[1,B]),2)/B;
    
    %Significant neurons. 
    sig = pval<0.05 & I>0;
    
    save('TemporalInfo.mat','I','sig');
    cd(currdir); 
end

function I = tempInfoOneNeuron(ratebylap)
    [nLaps,nBins,~] = size(ratebylap);
    
    p = 1/nBins;                                    %Probability occupancy in time (uniform).
    lambda_i = mean(4*ratebylap);                   %Mean rate at this time bin. (multiply by 4 because 0.25 time bins)
    lambda = sum(sum(4*ratebylap))/(nLaps*nBins);   %Mean rate. 
    
    %Temporal information.
    I = nansum(p.*lambda_i.*log2(lambda_i./lambda));
end