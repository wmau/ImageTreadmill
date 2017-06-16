function DecodePlacewithTemporalLookup(md,varargin)
%
%
%

%%
    path = md.Location;
    cd(path); 
    
    p = inputParser; 
    p.addRequired('md');
    load(fullfile(path,'FinalOutput.mat'),'NumNeurons');
    p.addParameter('neurons',1:NumNeurons,@(x) isnumeric(x)); 
    
    p.parse(md,varargin{:});
    
    neurons = p.Results.neurons;
    B = 1000;
    purple = [.58 .44 .86];
    teal = [0 .5 .5];
    nNeurons = length(neurons);
    
%% 
    %Linearize position and bin it. Also extract trial timestamps.
    [~,~,~,~,X,parsed] = LinearizedPF_raster(md,'plotit',false,'neurons',[]);
    
    %Determine the number of trials and longest trial. 
    trialDurs = tabulate(parsed.trial); 
    nTrials = max(parsed.trial);
    longestTrialDuration = max(trialDurs(:,2));
    nSpatialBins = max(X); 
%%
    %Load fluorescence data. 
    load(fullfile(path,'Pos_align.mat'),'RawTrace'); 
    
    %Preallocate.
    position = nan(nTrials,longestTrialDuration); 
    trace = nan(nTrials,longestTrialDuration,NumNeurons);
    lookup = nan(NumNeurons,longestTrialDuration);
    
    %Get position and trace information. 
    for r=1:nTrials
        thisTrial = find(parsed.trial==r);
        dur = length(thisTrial);
        
        %Get position and throw it into the matrix. 
        position(r,1:dur) = X(thisTrial);
        
        for n=1:nNeurons
            trace(r,1:dur,neurons(n)) = RawTrace(neurons(n),thisTrial);
        end
    end
    
%% Make lookup table. 
    %Simply average across trials. 
    for n=1:nNeurons
        lookup(neurons(n),:) = nanmean(trace(:,:,neurons(n)));
    end
    
%% 
    [decoded,observed] = deal(nan(B,nSpatialBins,NumNeurons));
    pLapSample = 0.5;
    nLapSample = round(pLapSample*nTrials);
    for n=1:nNeurons
        for i=1:B
            laps = randsample(1:nTrials,nLapSample);
            for b=1:nSpatialBins
                [~,t] = find(position(laps,:)==b);
                FDuringOccTime = lookup(neurons(n),t);
                decoded(i,b,neurons(n)) = mean(FDuringOccTime);
            end
        end
    end
    keyboard; 
end