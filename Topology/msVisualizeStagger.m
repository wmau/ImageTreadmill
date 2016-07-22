function msVisualizeStagger(mapMD,md,neuron,A,varargin)
%
%
%

%%
    p = inputParser;
    p.addRequired('mapMD',@(x) isstruct(x));
    p.addRequired('md',@(x) isstruct(x));
    p.addRequired('neuron',@(x) isnumeric(x) && isscalar(x));
    p.addRequired('A',@(x) islogical(x) || isnumeric(x));
    p.addParameter('edgelist',find(A(:,neuron))',@(x) isnumeric(x));
    
    p.parse(mapMD,md,neuron,A,varargin{:});
    el = p.Results.edgelist;
    nTriggers = length(el);
    nSessions = length(md);
    
    iDir = pwd;
       
    DATA = CompileMultiSessionData(md,{'ratebylap','ft','ttl','t'});
    targets = msMatchCells(mapMD,md,neuron);
    triggers = msMatchCells(mapMD,md,el);
    
    %Line format for second neuron raster.
    lead.Color = 'r';
    lead.LineWidth = 1;
    
    %Line format for first neuron raster, all spikes, transparent. 
    lag.Color = [0 1 0 0.3];    %Transparent green.
    lag.LineWidth = 2;
    
    %Line format for first neuron raster, only spikes immediately
    %preceding those of neuron two. 
    immediatelag.Color = 'g';
    immediatelag.LineWidth = 2;
    
    nTicks = 6;
    i = 1;
     
    for e=1:nTriggers
        f = figure('Position',[80 240 380*nSessions 570]);
        
        goodSessions = find(targets > 0 & triggers(e,:) > 0);
        
        for s=goodSessions
            nLaps = sum(DATA.ttl{s}.complete);
            
            %% Raster - Trigger
            subplot(3,2*nSessions,s*2-1);
            imagesc([0:DATA.t{s}],...
                [1:nLaps],...
                DATA.ratebylap{s}(:,:,triggers(e,s))); 
            colormap gray; title(['\color{green}Trigger \color{black}ROI #',num2str(triggers(e,s))]);
            if s==1, ylabel('Laps'); end

            %% Raster - target
            subplot(3,2*nSessions,s*2);
            imagesc([0:DATA.t{s}],...
                [1:nLaps],...
                DATA.ratebylap{s}(:,:,targets(s)));
            colormap gray; title(['\color{red}Target \color{black}ROI #',num2str(targets(s))]);

            %% Tick raster
            inds = DATA.ttl{s}.inds(logical(DATA.ttl{s}.complete));
            inds(:,2) = inds(:,1) + 20*DATA.t{s}-1;
            triggerRaster = buildRaster(inds,DATA.ft{s},triggers(e,s));
            targetRaster = buildRaster(inds,DATA.ft{s},targets(s));

            [immediateRaster,d] = stripRaster(triggerRaster,targetRaster);

            subplot(3,2*nSessions,s*2+nSessions*2-1:s*2+nSessions*2);
            plotSpikeRaster(triggerRaster,'PlotType','vertline',...
                'LineFormat',lag,'TimePerBin',0.05,'SpikeDuration',0.05); 
            hold on;
            plotSpikeRaster(immediateRaster,'PlotType','vertline',...
                'LineFormat',immediatelag,'TimePerBin',0.05,'SpikeDuration',0.05); 
            plotSpikeRaster(targetRaster,'PlotType','vertline',...
                'LineFormat',lead,'TimePerBin',0.05,'SpikeDuration',0.05);   
                ax = gca; 
                ax.Color = 'k';
                ax.XTick = linspace(ax.XLim(1),ax.XLim(2),nTicks);
                ax.XTickLabel = linspace(0,DATA.t{s},nTicks);
                ax.YTick = [1:5:nLaps];
                set(gca,'ticklength',[0 0]);
                hold off; if s==1, ylabel('Laps'); end

            %% Target to trigger/treadmill latency
            bothActiveLaps = find(any(immediateRaster,2));
            TMAlignedOnsets = [];
            for l=bothActiveLaps'
                %Get the onset times of each neuron. 
                TMAlignedOnsets = [TMAlignedOnsets find(targetRaster(l,:))];    
            end
            %[~,TMalignedOnsets] = find(targetRaster);

            %Divide by frame rate. 
            TMAlignedOnsets = TMAlignedOnsets./20;

            %Spread of responses relative to treadmill start. 
            treadmillOffsetSpread = mad(TMAlignedOnsets,1);

            subplot(3,2*nSessions,s*2+nSessions*4-1:s*2+nSessions*4);
            histogram(TMAlignedOnsets,[0:0.25:10],'normalization','probability',...
                'facecolor','k');
            hold on;      

            %Get spread. 
            cellOffsetSpread = mad(d,1);

            %Ratio between cell-to-cell vs cell-to-treadmill.
            ratio =  cellOffsetSpread / treadmillOffsetSpread;

            %Histogram.
            histogram(-d,[0:0.25:10],'normalization','probability',...
                'facecolor','y');
            title(['TT Score = ',num2str(ratio)]);
            legend({'Treadmill','Trigger'});
            xlabel('Latency from Target [s]'); 
            if s==1, ylabel('Proportion of Spikes'); end
            set(gca,'linewidth',1.5);
            hold off;
        end
            
        set(f,'PaperOrientation','landscape');
        set(f,'PaperUnits','normalized');
        set(f,'PaperPosition',[0 0 1 1]);
    end
        
end