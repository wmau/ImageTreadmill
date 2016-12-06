function spatInfo(md)
% spatInfo(md)
%
%   Calculates the Shannon mutual information I(X,K) between the random
%   variables spike count [0,1] and position via the equations: 
%
%   (1) I_pos(xi) = sum[k>=0](P_k|xi * log2(P_k|xi / P_k)) 
%
%   (2) MI = sum[i=1->N](P_xi * I_pos(xi)
%
%   where
%       P_xi is the probability the mouse is in pixel xi,
%       RunOccMap./sum(RunOccMap(:)
%       
%       P_k is the probability of observe k spikes,
%       sum(FT(neuron,:),2)/size(FT,2)
%
%       P_k|xi is the conditional probability of observing k spikes in
%       pixel xi, TMap_unsmoothed
%

%% Load.
    cd(md.Location);
    load('PlaceMaps.mat','TMap_unsmoothed','RunOccMap','isrunning','cmperbin',...
        'frames_use_ind');
    load('Pos_align.mat','FT');
    
%% Set up variables. 
    %Get good frames and good pixels. 
    okframes = frames_use_ind & isrunning;      %Running, but not on treadmill.
        
    %Number of frames and neurons. 
    nFrames = sum(okframes); 
    nNeurons = size(FT,1);
    
    %Get dwell map. 
    P_x = RunOccMap(:)./nFrames;
    okpix = RunOccMap(:) > 4 & P_x < .05;       %Dwell must be for more than 4 frames, but less than 5% of total session.
    P_x = P_x(okpix);                           %Only take good pixels.
    
    %Get probability of spiking and not spiking for each neuron.
    P_k1 = sum(FT(:,okframes),2)./nFrames;
    P_k0 = 1-P_k1;
    
%% Compute information metrics. 
    %Preallocate.
    P_1x = cell(1,nNeurons);        %Probability of spike given location.
    P_0x = cell(1,nNeurons);        %Probability of ~spike given location.
    Ipos = cell(1,nNeurons);        %Positional information vector for each neuron.
    MI = zeros(1,nNeurons);         %Mutual information for each neuron.
    for n=1:nNeurons
        %Get probability of spike given location, TMap, only taking good
        %pixels. 
        P_1x{n} = TMap_unsmoothed{n}(:);    P_1x{n} = P_1x{n}(okpix);
        P_0x{n} = 1-TMap_unsmoothed{n}(:);  P_0x{n} = P_0x{n}(okpix);
        
        %Compute positional information for k=1 and k=0.
        I_k1 = P_1x{n}.*log(P_1x{n}./P_k1(n));
        I_k0 = P_0x{n}.*log(P_0x{n}./P_k0(n)); 
        
        %Sum these to make true positional information. 
        Ipos{n} = I_k1 + I_k0;
        
        %Compute mutual information.
        MI(n) = nansum(P_x.*Ipos{n});   %bits. 
        
        %Compute information content per sec or spk.
        Isec(n) = nansum(P_1x{n}.*P_x.*log2(P_1x{n}./P_k1(n))); %bits/sec
        Ispk(n) = Isec(n) ./ P_k1(n);                           %bits/spk    
    end
    
    save('SpatialInfo.mat','MI','Isec','Ispk','Ipos','okframes','okpix');
end