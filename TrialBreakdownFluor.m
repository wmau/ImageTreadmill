function TrialBreakdownFluor(md,neurons,nBlocks)
%
%
%

%%
    path = md.Location;
    cd(path);
    
    load('TimeCells.mat','T');
    load('TreadmillTraces.mat','DFDTTrdmll');
    [nTrials,nFrames,~] = size(DFDTTrdmll); 
    blockLims = floor(linspace(0,nTrials,nBlocks+1));
   
    keepgoing = true;
    thisNeuron = 1;
    t = linspace(0,T,nFrames);
    while keepgoing
        f = figure(neurons(thisNeuron));
        f.Position = [580 100 750 870];
        
        for b=1:nBlocks
            span = blockLims(b)+1:blockLims(b+1);
            
            nSpan = length(span);
            %colors = parula(nSpan); 
            colors = zeros(nSpan,3);
            
            subplots(b) = subplot_auto(nBlocks,b); hold on;
            for s=1:nSpan
                plot(t,DFDTTrdmll(span(s),:,neurons(thisNeuron))','color',colors(s,:));
            end
           
            plot(t,mean(DFDTTrdmll(span,:,neurons(thisNeuron))),'color','k',...
                'linewidth',3);
            
            title(['Trials ',num2str(span(1)),' - ',num2str(span(end))]);  
            set(gca,'tickdir','out','linewidth',4,'fontsize',12);
        end
        
        yLims = [min([subplots.YLim]) max([subplots.YLim])];
        set(subplots,'ylim',yLims);
        
        [keepgoing,thisNeuron] = scroll(thisNeuron,length(neurons),f);
        close all;
    end
end