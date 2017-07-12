function TrialBreakdownFluor_time(md,neurons,nBlocks)
%TrialBreakdownFluor_time(md,neurons,nBlocks)
%
%   Splits trials into n blocks chronologically, then plots fluorescence
%   traces by trial. 
%
%   INPUTS
%       md: session entry.
%
%       neurons: vector of neurons to scroll through. 
%
%       nBlocks: scalar indicating number of trial blocks.
%

%% Set up.
    path = md.Location;
    cd(path);
    
    load('TimeCells.mat','T');
    load('TreadmillTraces.mat','DFDTTrdmll');
    [nTrials,nFrames,~] = size(DFDTTrdmll); 
    blockLims = floor(linspace(0,nTrials,nBlocks+1));
   
%% Loop through neurons. 
    keepgoing = true;
    thisNeuron = 1;
    t = linspace(0,T,nFrames);
    while keepgoing
        f = figure(neurons(thisNeuron));
        f.Position = [580 100 750 870];
        
        for b=1:nBlocks
            span = blockLims(b)+1:blockLims(b+1);   %Trial block limits. 
            
            %Number of trials per block. 
            nSpan = length(span);
            %colors = parula(nSpan); 
            colors = zeros(nSpan,3);
            
            %Plot individual trials. 
            subplots(b) = subplot_auto(nBlocks,b); hold on;
            for s=1:nSpan
                plot(t,DFDTTrdmll(span(s),:,neurons(thisNeuron))','color',colors(s,:));
            end
           
            %Plot mean of that trial block. 
            plot(t,mean(DFDTTrdmll(span,:,neurons(thisNeuron))),'color',[0 .5 .5],...
                'linewidth',3);
            
            %Label. 
            title(['Trials ',num2str(span(1)),' - ',num2str(span(end))]);  
            set(gca,'tickdir','out','linewidth',4,'fontsize',12);
            axis tight;
        end
        
        %Normalize axis limits. 
        yLims = [min([subplots.YLim]) max([subplots.YLim])];
        set(subplots,'ylim',yLims);
        
        [keepgoing,thisNeuron] = scroll(thisNeuron,length(neurons),f);
        close all;
    end
end