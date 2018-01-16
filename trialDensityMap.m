function [TCtrialDensity,neurons,order] = trialDensityMap(md,varargin)
%
%
%

%%
    neurons = getTimeCells(md);
    load(fullfile(md.Location,'TimeCells.mat'),'TodayTreadmillLog','ratebylap',...
        'curves'); 
    complete = logical(TodayTreadmillLog.complete);
    nTrialBlocks = sum(complete); 
    
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('nTrialBlocks',nTrialBlocks,@(x) isscalar(x)); 
    p.addParameter('neurons',neurons,@(x) isnumeric(x));
    p.addParameter('order',[]); 
    
    p.parse(md,varargin{:});
    nTrialBlocks = p.Results.nTrialBlocks; 
    neurons = p.Results.neurons;
    order = p.Results.order;
    
%%
    ratebylap = ratebylap(complete,:,:);
    [nRuns,~,nNeurons] = size(ratebylap); 
    ratebylap = ratebylap > 0;
    sig = curves.sig; 
    
    trialBlockLims = linspace(0,nRuns,nTrialBlocks+1); 
    trialBlockLims = floor(trialBlockLims); 
    
%%
    trialDensity = nan(nNeurons,nTrialBlocks); 
    for neuron=1:nNeurons
        goodbins = sig{neuron};
        for b=1:nTrialBlocks
            block = trialBlockLims(b)+1:trialBlockLims(b+1);
            
            trialDensity(neuron,b) = mean(any(ratebylap(block,goodbins,neuron),2));
        end
    end
    
    TCtrialDensity = trialDensity(neurons,:);
    TCtrialDensity = TCtrialDensity./max(TCtrialDensity,[],2);
    neverActive = ~any(TCtrialDensity,2); 
    [~,ind] = max(TCtrialDensity,[],2);
    TCtrialDensity(neverActive,:) = nan;
    
    if isempty(order)
        skewness = getAllSkewnesses(md);
        [~,order] = sort(skewness(neurons));
    end
    
    %figure('position',[825     10   260   990]);
    h = imagesc(TCtrialDensity(order,:));
    set(h,'AlphaData',~isnan(TCtrialDensity(order,:)));
    %axis equal; 
    axis tight;
    xlabel('Trial blocks');
    ylabel('Cell #');
    set(gca,'fontsize',14,'xtick',[]);
    
end