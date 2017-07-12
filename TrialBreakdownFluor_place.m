function TrialBreakdownFluor_place(md,neurons,nBlocks,varargin)
%TrialBreakdownFluor_place(md,neurons,nBlocks,varargin)
%
%   Splits trials into n blocks chronologically, then plots fluroescence
%   traces by trial. 
%   
%   INPUTS
%       md: session entry.
%
%       neurons: vector of neurons to scroll through. 
%
%       nBlocks: scalar indicating number of trial blocks. 
%
%       NAME,VALUE
%           noTreadmill: logical indicating whether or not to include the
%           treadmill in the plot. 
%

%% Set up.
    p = inputParser;
    p.addRequired('md',@(x) isstruct(x));
    p.addRequired('neurons',@(x) isnumeric(x));
    p.addRequired('nBlocks',@(x) isnumeric(x)); 
    p.addParameter('noTreadmill',true,@(x) islogical(x));
    
    p.parse(md,neurons,nBlocks,varargin{:});
    
    noTreadmill = p.Results.noTreadmill;
    
    path = md.Location;
    cd(path); 
    nBins = 80;
    
    load('Pos_align.mat','DFDTtrace'); 

%% Linearize space. 
    [~,~,~,~,binOcc,parsed] = LinearizedPF_raster(md,'plotit',false,...
        'neurons',1,'nBins',nBins);
    
    nTrials = max(parsed.trial); 
    blockLims = floor(linspace(0,nTrials,nBlocks+1)); 
    
%% Loop through neurons. 
    keepgoing = true;
    thisNeuron = 1;
    x = 1:nBins;  
    while keepgoing
        f = figure(neurons(thisNeuron));
        f.Position = [580 100 750 870];
        
        for b=1:nBlocks
            span = blockLims(b)+1:blockLims(b+1);
            
            %Number of trials per block. 
            nSpan = length(span);
            F = nan(nSpan,nBins);
            subplots(b) = subplot_auto(nBlocks,b); hold on;
            for s=1:nSpan
                binOccTrial = binOcc(parsed.trial==span(s)); 
                traceTrial = DFDTtrace(neurons(thisNeuron),parsed.trial==span(s)); 
                
                %Build vector, taking the mean of fluorescence while in
                %that bin. 
                F(s,:) = accumarray(binOccTrial',traceTrial',[nBins,1],@mean);
                plot(x,F,'color','k');       
            end
            
            %Plot mean of that trial block. 
            plot(x,mean(F),'color','k','linewidth',3);
            title(['Trials ',num2str(span(1)),' - ',num2str(span(end))]); 
            set(gca,'tickdir','out','linewidth',4,'fontsize',12);
            
            %If no treadmill, remove the treadmill bins. 
            if noTreadmill
                xlim([22 80]);
                set(gca,'xtick',[22 80],'xticklabel',[0 140]); 
            end
        end
        
        %Normalize axis limits.  
        yLims = [min([subplots.YLim]) max([subplots.YLim])];
        set(subplots,'ylim',yLims);
        
        [keepgoing,thisNeuron] = scroll(thisNeuron,length(neurons),f);
        close all;
    end
end