function trialDensityMap(md,varargin)
%
%
%

%%
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x)); 
    p.addParameter('nTrialBlocks',10,@(x) isscalar(x)); 
    
    p.parse(md,varargin{:});
    nTrialBlocks = p.Results.nTrialBlocks; 
    
%%
    load(fullfile(md.Location,'TimeCells.mat'),'TodayTreadmillLog','ratebylap',...
        'curves'); 
    complete = logical(TodayTreadmillLog.complete);
    ratebylap = ratebylap(complete,:,:);
    [nRuns,~,nNeurons] = size(ratebylap); 
    sig = curves.sig; 
    
    trialBlockLims = linspace(0,nRuns,nTrialBlocks+1); 
    trialBlockLims = floor(trialBlockLims); 
    
%%
    trialDensity = nan(nNeurons,nTrialBlocks); 
    for neuron=1:nNeurons
        goodbins = sig{neuron};
        for b=1:nTrialBlocks
            block = trialBlockLims(b)+1:trialBlockLims(b+1);
            
            trialDensity(neuron,b) = mean(sum(ratebylap(block,goodbins,neuron),2));
        end
    end
    
    TimeCells = getTimeCells(md);
    TCtrialDensity = trialDensity(TimeCells,:);
    TCtrialDensity = TCtrialDensity./max(TCtrialDensity,[],2);
    [~,ind] = max(TCtrialDensity,[],2);
    [~,order] = sort(ind);
    
    figure('position',[825     10   260   990]);
    imagesc(TCtrialDensity(order,:));
    axis equal; axis tight;
    xlabel('Trial blocks');
    ylabel('Cell #');
    set(gca,'fontsize',14,'xtick',[]);
    
end